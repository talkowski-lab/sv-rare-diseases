version 1.0

workflow Episignature {
    input {
        Array[File] combined_beds
        File? signature_bed
        File? case_signature
        File? control_signature
        String? condition_name
        String prefix
 
    }
    meta {
        description: "Episignature analysis workflow, the input is a list of combined bed files"
    }

    if (defined(signature_bed) && defined(case_signature) && defined(control_signature) && defined(condition_name)) {
        call CustomEpisignaturePreProcessBeds {
            input:
                signature_bed = select_first([signature_bed]),
                combined_beds = select_all(combined_beds)
        }
        call CustomEpisignatureMatching {
            input:
                episignature_bed_files = CustomEpisignaturePreProcessBeds.episignature_bed_files,
                signature_bed = select_first([signature_bed]),
                case_signature = select_first([case_signature]),
                control_signature = select_first([control_signature]),
                condition_name = select_first([condition_name]),
                output_prefix = prefix
        }
    }
    if (! (defined(signature_bed) && defined(case_signature) && defined(control_signature) && defined(condition_name))) {

        call PreProcessBeds {
            input:
                combined_beds = select_all(combined_beds)
        }

        call EpisignatureMatching {
            input:
                episignature_bed_files = PreProcessBeds.episignature_bed_files,
                output_prefix = prefix
        }
    }

    output {
        Array[File] episignature_bed_files = select_first([CustomEpisignaturePreProcessBeds.episignature_bed_files ,PreProcessBeds.episignature_bed_files])
        #File? assignment = EpisignatureMatching.assignment
        File results = select_first([CustomEpisignatureMatching.results,EpisignatureMatching.results])
    }
}


task PreProcessBeds {
    input {
        Array[File] combined_beds
    }


    command <<<
        set -euxo pipefail
        mkdir /episignature/NSBEpi/hg38_episignature_cordinates_with_chr/

        for file in /episignature/NSBEpi/hg38_episignature_cordinates/*bed; \
            do base=$(basename $file); \
            sed 's/^/chr/' $file > /episignature/NSBEpi/hg38_episignature_cordinates_with_chr/$base ; \
        done

        # Creating the output folder if it doesn't exist
        mkdir -p combined_beds/
        mkdir -p combined_beds/episignature_loci_v2/
        
        for file in ~{sep=" " combined_beds}; do
            file2_basename=$(basename "$file" .combined.bed.gz)
            echo "Processing: $file2_basename"
            
            # Decompress and sort the input file
            gunzip -c "$file" | sort -k1,1 -k2,2n > "combined_beds/${file2_basename}.sorted.bed"
            
            # Intersect with each reference file
            for file1 in /episignature/NSBEpi/hg38_episignature_cordinates_with_chr/*.bed; do
                if [ ! -f "$file1" ]; then
                    echo "Warning: Reference file not found: $file1"
                    continue
                fi
                
                file1_basename=$(basename "$file1" .bed)
                output_file="combined_beds/episignature_loci_v2/${file1_basename}_${file2_basename}.bed"
                
                bedtools intersect -sorted -a "$file1" -b "combined_beds/${file2_basename}.sorted.bed" \
                    -loj -wa -wb > "$output_file"
            done
            
            # Clean up temporary file
            rm "combined_beds/${file2_basename}.sorted.bed"
        done

    >>>

    output {
        Array[File] episignature_bed_files = glob("combined_beds/episignature_loci_v2/*.bed")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/episignatures:v0"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}

task CustomEpisignaturePreProcessBeds {
    input {
        Array[File] combined_beds
        String condition_name
        File signature_bed
    }


    command <<<
        set -euxo pipefail

        # Creating the output folder if it doesn't exist
        mkdir -p combined_beds/
        mkdir -p combined_beds/episignature_loci_v2/
        
        for file in ~{sep=" " combined_beds}; do
            file2_basename=$(basename "$file" .combined.bed.gz)
            echo "Processing: $file2_basename"
            
            # Decompress and sort the input file
            gunzip -c "$file" | sort -k1,1 -k2,2n > "combined_beds/${file2_basename}.sorted.bed"
            
            # Intersect with each reference file                
            output_file="combined_beds/episignature_loci_v2/~{condition_name}_${file2_basename}.bed"
            
            bedtools intersect -sorted -a ~{signature_bed} -b "combined_beds/${file2_basename}.sorted.bed" \
                -loj -wa -wb > "$output_file"
            
            # Clean up temporary file
            rm "combined_beds/${file2_basename}.sorted.bed"
        done

    >>>

    output {
        Array[File] episignature_bed_files = glob("combined_beds/episignature_loci_v2/*.bed")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/episignatures:v0"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}


task EpisignatureMatching {
    input {
        Array[File] episignature_bed_files
        String output_prefix

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        python - --bedfiles "~{sep=' ' episignature_bed_files}" \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import argparse
        import os
        import pandas as pd
        import pickle
        import numpy as np
        from sklearn import svm
        from sklearn.metrics import accuracy_score
        from scipy import stats
        import json
        from sklearn.model_selection import cross_validate, StratifiedKFold
        from sklearn.metrics import (accuracy_score, precision_score, recall_score, 
                             f1_score, roc_auc_score, confusion_matrix)
        from sklearn.metrics import make_scorer, f1_score, precision_score

        def read_bed_files(sample_bed_files, epi_signatures_all):
            bed_files = {}
            sorted_files = sorted(sample_bed_files)
            n_samples = int(len(sorted_files)/len(epi_signatures_all))
            sample_names = sorted_files[0:n_samples]
            sample_names = [s.split("/")[-1] for s in sample_names]
            adjusted_sample_names = [s.split('.')[0] for s in sample_names] # drop the suffix
            adjusted_sample_names = ['_'.join(s.split('_')[1:]) for s in adjusted_sample_names] # drop the disorder_ prefix
            # Iterate over each file and read the bedfiles 
            for file_path in sorted_files:
                if file_path.endswith('.bed'):
                    print(file_path)
                    # Read the contents of the bedmethyl file using pandas (first 3 columsn chr start end and methylation column 13th position)
                    bed_df = pd.read_csv(file_path, sep='\t', header=None)
                    bed_df = (bed_df[[0, 1, 2, 11]]).reset_index(drop=True)
                    bed_df = bed_df.rename(columns={0: 'Chr', 1: 'Start', 2: 'End', 11 : 'Methylation'})
                    file_basename = os.path.splitext(file_path)[0]  
                    # correct the formatting of methylation values
                    bed_df['Methylation'] = bed_df['Methylation'].replace('.', np.nan)
                    bed_df['Methylation'] = pd.to_numeric(bed_df['Methylation'] , errors='coerce')
                    bed_df['Methylation'] = bed_df['Methylation']/100
                    #sample_order.append(file_basename)
                    bed_files[file_basename] = bed_df
            return bed_files, n_samples, adjusted_sample_names


        def run_SVM(epi_signatures_all_with_sample, n_of_samples):
            filter_df = pd.DataFrame()
            pred_df = pd.DataFrame()
            order_disorders = []

            for disorder in epi_signatures_all_with_sample:
                if (disorder != 'Controls'):

                    X = epi_signatures_all_with_sample[disorder].T
                    y = X.index.tolist()

                    if (disorder != 'MRXCJS'):
                        new_y = [1 if item == str(disorder) else 0 for item in y]
                    else:
                        new_y = [1 if item == 'MRXSCJ' else 0 for item in y]

                    # Convert 'NA' values to NaN and Fill NaN values with column averages
                    X = X.replace('NA', np.nan)
                    X = X.astype(float)
                    X = X.fillna(X.mean())


                    #take for training always the first 35 samples representing the illumina data
                    #take for testing always our nanopore samples
                    X_train = X.iloc[:-n_of_samples, :]
                    X_test = X.iloc[-n_of_samples:, :]

                    y_train = new_y[:-n_of_samples]
                    y_test = new_y[-n_of_samples:]

                    #weights for classes
                    #class_weight='balanced'
                    class_weights = {0: 1, 1: 10} 

                    # Create SVM classifiers 
                    linear_classifier = svm.SVC(kernel='linear', class_weight=class_weights, probability=True)

                    if (X_train.shape[0] > 0) and (len(np.unique(y_train))>1):

                        #Train the classifiers for the current disorder
                        linear_classifier.fit(X_train, y_train)

                        # Make predictions using the trained classifiers
                        linear_pred = linear_classifier.predict(X_test)

                        #store the prediction results of every SVM
                        new_row_df = pd.DataFrame([linear_pred])
                        pred_df =  pd.concat([pred_df,new_row_df], ignore_index=True)
                        #store the order in which the disorder specific SVMs were trained and tested
                        order_disorders.append(disorder)

                        #store decision function value for every SVM - this will be used to determine the samples correct class
                        decision_values = linear_classifier.decision_function(X_test)
                        new_row_df = pd.DataFrame([decision_values])
                        filter_df = pd.concat([filter_df,new_row_df], ignore_index=True)
            return filter_df, pred_df, order_disorders

        def find_assignments(filter_df, samples, order_disorders):
            n_of_samples=len(set(samples))
            assignments = []
            for sample in range(0, n_of_samples):
                assigned_values = []
                assigned_disorders = []
                i=0
                max_value = filter_df.iloc[:,sample].max()
                row_name = filter_df.iloc[:,sample].idxmax()

                for svm_value in filter_df.iloc[:,sample]:
                    if svm_value >= 0.1:
                        assigned_values.append(svm_value)
                        assigned_disorders.append(order_disorders[i])
                    i+=1

                # Print the maximum value and row name
                # print(samples[sample],end="\t")
                if (max_value == -1000 or max_value < 0.1):
                    assigned_disorders = ['Control']

                if len(assigned_values) == 1:
                    assignments.append([samples[sample], str(round(max_value,3)), assigned_disorders[0]])
                elif len(assigned_values) >1:
                    assignments.append([samples[sample], ','.join([str(round(x,3)) for x in assigned_values]),','.join(assigned_disorders)])
                else:
                    assignments.append([samples[sample], str(round(max_value,3))])
            return assignments

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--bedfiles',
                                type=str)

            parser.add_argument('--output_file',
                                type=str)


            args = parser.parse_args()
            # path to pickle file containing a dictionary with illumina data preprocessed
            file_path = '/episignature/NSBEpi/data/no_strand_all_points_dict.pickle'
            with open(file_path, 'rb') as file:
                epi_signatures_all = pickle.load(file)
            # load sample bed files
            input_path = args.bedfiles.split(" ")
            result, n_of_samples, samples = read_bed_files(input_path, epi_signatures_all)
            print(samples)

            #link nanopore methylation data to the correct disorder
            epi_signatures_all_with_sample = epi_signatures_all.copy()
            epi_signatures_all_with_sample['MRXCJS'] = epi_signatures_all_with_sample['MRXSCJ']
            del epi_signatures_all_with_sample['MRXSCJ']
            for disorder in epi_signatures_all_with_sample:
                epi_signatures_all_with_sample[disorder] = epi_signatures_all_with_sample[disorder].iloc[:, 2:-1]
                for file_basename in result:
                    if file_basename.startswith(disorder):
                        epi_signatures_all_with_sample[disorder] = epi_signatures_all_with_sample[disorder].assign(**{file_basename: list(result[file_basename]['Methylation'])})
            # filter and predict
            filter_df, pred_df, order_disorders = run_SVM(epi_signatures_all_with_sample, n_of_samples )
            # find assignments
            assignments = find_assignments(filter_df, samples, order_disorders)
            # write output file
            with open(args.output_file + ".assignment.txt", "w") as out_f:
                out_f.write("Sample\tMax_value\tAssigned_disorder\n")
                for ass in assignments:
                    out_f.write("\t".join(ass) + "\n")

            #determine for every sample the disorder with the highest decision function value
            filter_df.columns = samples[:n_of_samples]
            filter_df.index = order_disorders
            filter_df.to_csv(args.output_file + ".episignature.results.tsv", sep="\t")

        if __name__ == "__main__":
            main()
        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/episignatures:v0"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }

    output {
        File assignment = "~{output_prefix}.assignment.txt"
        File results = "~{output_prefix}.episignature.results.tsv"
    }
}


task CustomEpisignatureMatching {
    input {
        Array[File] episignature_bed_files
        
        File signature_bed
        File case_signature
        File control_signature
        String output_prefix
        String condition_name

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        mkdir -p unzipped_beds

        unzipped_files=()

        # if files are gzipped, unzip them first
        for bed in ~{sep=" " episignature_bed_files}; do
            bed_basename=$(basename "$bed" .gz)
            unzipped_path="unzipped_beds/$bed_basename"
            if [[ "$bed" == *.gz ]]; then
                gunzip -c "$bed" > "$unzipped_path"
            else
                cp "$bed" "$unzipped_path"
            fi
            unzipped_files+=("$unzipped_path")
        done

        python - --bedfiles "${unzipped_files[*]}" \
                 --case_signature ~{case_signature} \
                 --control_signature ~{control_signature} \
                 --condition_name ~{condition_name} \
                 --signature_bed ~{signature_bed} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import argparse
        import os
        import pandas as pd
        import pickle
        import numpy as np
        from sklearn import svm
        from sklearn.metrics import accuracy_score
        from scipy import stats
        from sklearn.model_selection import cross_validate, StratifiedKFold
        from sklearn.metrics import (accuracy_score, precision_score, recall_score, 
                             f1_score, roc_auc_score, confusion_matrix)
        from sklearn.metrics import make_scorer, f1_score, precision_score
        import json

        def read_bed_files(sample_bed_files):
            bed_files = {}
            sorted_files = sorted(sample_bed_files)
            sample_names = [s.split("/")[-1] for s in sorted_files]
            adjusted_sample_names = [s.split('.')[0] for s in sample_names] # drop the suffix
            adjusted_sample_names = ['_'.join(s.split('_')[1:]) for s in adjusted_sample_names] # drop the disorder_ prefix
            # Iterate over each file and read the bedfiles 
            for i, file_path in enumerate(sorted_files):
                if file_path.endswith('.bed'):
                    # Read the contents of the bedmethyl file using pandas (first 3 columsn chr start end and methylation column 13th position)
                    bed_df = pd.read_csv(file_path, sep='\t', header=None)
                    bed_df = (bed_df[[3, 12]]).reset_index(drop=True)
                    bed_df = bed_df.rename(columns={3: 'ID_REF', 12 : 'Methylation'})
                    file_basename = adjusted_sample_names[i]
                    # correct the formatting of methylation values
                    bed_df['Methylation'] = bed_df['Methylation'].replace('.', np.nan)
                    bed_df['Methylation'] = pd.to_numeric(bed_df['Methylation'] , errors='coerce')
                    bed_df['Methylation'] = bed_df['Methylation']/100
                    bed_df = bed_df.set_index("ID_REF")
                    bed_files[file_basename] = bed_df
            return bed_files, adjusted_sample_names

        def evaluate_svm_with_cv(X_train, y_train, class_weights, cv_folds):
            """
            Evaluate SVM model performance using cross-validation.
            """
            # Create classifier with specified class weights
            classifier = svm.SVC(kernel='linear', class_weight=class_weights, probability=True)
            
            # Define scoring metrics
            scoring = {
                'accuracy': 'accuracy',
                'precision': make_scorer(precision_score, pos_label=1, zero_division=0),
                'recall': make_scorer(recall_score, pos_label=1, zero_division=0),
                'f1': make_scorer(f1_score, pos_label=1, zero_division=0),
                'roc_auc': 'roc_auc'
            }
            
            # Use StratifiedKFold to maintain class distribution
            cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
            
            # Perform cross-validation
            cv_results = cross_validate(
                classifier,
                X_train,
                y_train,
                cv=cv,
                scoring=scoring,
                return_train_score=True,
                n_jobs=-1
            )
            cv_results_new = {}
            for key, value in cv_results.items():
                cv_results_new[key] = value.tolist()
                        
            return cv_results_new

        def run_SVM_single(df_training, df_testing, cv_folds = 5):            
            REF_ID_list = df_training.index.tolist()
            df_testing = df_testing.reindex(REF_ID_list)
            df_train_test = pd.concat([df_training, df_testing], axis = 1)
            assert df_train_test.shape[0] == df_testing.shape[0]
            X_train_test = df_train_test.T
            # Convert 'NA' values to NaN and Fill NaN values with column averages (Warning: may revisit)
            X_train_test = X_train_test.replace('NA', np.nan)
            X_train_test = X_train_test.astype(float)
            X_train_test = X_train_test.fillna(X_train_test.mean())
            
            # split train and test
            X_train = X_train_test.loc[df_training.columns, :]
            y_train_index = X_train.index.tolist()
            y_train = [1 if item.startswith("Case") else 0 for item in y_train_index]


            X_test = X_train_test.loc[df_testing.columns, :]


            #weights for classes
            #class_weights = {0: 1, 1: 10} 
            class_weights = "balanced"

            # Create SVM classifiers 
            linear_classifier = svm.SVC(kernel='linear', class_weight=class_weights, probability=True)

            # cv
            cv_results = evaluate_svm_with_cv(X_train, y_train, class_weights, cv_folds)

            #Train the classifiers for the current disorder
            linear_classifier.fit(X_train, y_train)

            # Make predictions using the trained classifiers
            linear_pred = linear_classifier.predict(X_test)

            #store the prediction results of every SVM
            pred_df = pd.DataFrame([linear_pred])
            
            #store decision function value for every SVM - this will be used to determine the samples correct class
            decision_values = linear_classifier.decision_function(X_test)
            filter_df = pd.DataFrame([decision_values])
                
            return filter_df, pred_df, cv_results

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--bedfiles',
                                type=str)
            parser.add_argument('--case_signature',
                                type=str)
            parser.add_argument('--control_signature',
                                type=str)
            parser.add_argument('--condition_name',
                                type=str)
            parser.add_argument('--signature_bed',
                                type=str)
            parser.add_argument('--output_file',
                                type=str)


            args = parser.parse_args()
            # load signature probelist
            Signature_prob = []
            with open(args.signature_bed, 'r') as fp:
                data = fp.readlines()
                for line in data:
                    itemlist = line[:-1].split('\t')
                    Signature_prob.append(itemlist[-1])
            # load training data
            case_prob_values = pd.read_csv(args.case_signature, sep = "\t", header = 0)
            case_data = case_prob_values[case_prob_values["ID_REF"].isin(Signature_prob)]
            case_data.set_index("ID_REF", inplace=True)
            case_data.columns = ["Case_" + c for c in case_data.columns]

            control_prob_values = pd.read_csv(args.control_signature, sep = "\t", header = 0)
            control_data = control_prob_values[control_prob_values["ID_REF"].isin(Signature_prob)]
            control_data.set_index("ID_REF", inplace=True)
            control_data.columns = ["Control_" + c for c in control_data.columns]

            # merge
            df_training = pd.concat([case_data, control_data], axis = 1)

            # load sample bed files
            input_path = args.bedfiles.split(" ")
            df_testing_dict, adjusted_sample_names = read_bed_files(input_path)
            print(adjusted_sample_names, df_testing_dict.keys())

            # filter and predict
            Results = []
            for sample_name, df_testing in df_testing_dict.items():
                filter_df, pred_df, cv_results = run_SVM_single(df_training, df_testing)
                Results.append([sample_name, filter_df.values[0][0], pred_df.values[0][0], args.condition_name, json.dumps(cv_results)])

            # write output file
            with open(args.output_file + ".prediction.txt", "w") as out_f:
                out_f.write("Sample\tDecision_value\tPrediction\tDisorder\tCV_results\n")
                for ass in Results:
                    out_f.write("\t".join([str(a) for a in ass]) + "\n")

        if __name__ == "__main__":
            main()
        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/episignatures:v0"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }

    output {
        File results = "~{output_prefix}.prediction.txt"

    }
}