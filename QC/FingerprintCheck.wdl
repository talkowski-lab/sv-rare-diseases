version 1.0

# Given two bam/cram files, use somalier to extract sites and compare their fingerprints

workflow FingerprintCheck {
    input {
        File bam1
        File bai1
        String prefix1="f1_"
        File bam2
        File bai2
        String prefix2="f2_"
        File reference
        File sites
    }

  call ExtractFingerprint as ExtractFingerprint1 {
    input:
      bam = bam1,
      bai = bai1,
      prefix = prefix1,
      sites = sites,
      reference = reference
  }

  call ExtractFingerprint as ExtractFingerprint2 {
    input:
      bam = bam2,
      bai = bai2,
      prefix = prefix2,
      sites = sites,
      reference = reference
  }

  call CheckFingerprint {
    input:
      fingerprint1 = ExtractFingerprint1.fingerprint,
      fingerprint2 = ExtractFingerprint2.fingerprint,
  }

  output {
    String relatedness = CheckFingerprint.relatedness
    File fingerprint1 = ExtractFingerprint1.fingerprint
    File fingerprint2 = ExtractFingerprint2.fingerprint
  }
}

task ExtractFingerprint {
    input {
        File bam
        File bai
        String prefix
        File reference
        File sites
        
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
        disk_gb: ceil(size(bam, "GB")+size(reference, "GB"))*5 + 20,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "brentp/somalier:v0.2.19"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    String base1 = basename(bam, ".cram")
    String base = basename(base1, ".bam") # remove .bam or .cram, whichever is present

    command <<<
        set -euxo pipefail
        
        mkdir extracted

        somalier extract -d extracted/ --sites ~{sites} -f ~{reference} ~{bam} --sample-prefix ~{prefix}

        # get generated filename
        filename=$(ls extracted/)
        
        mv extracted/${filename} ~{prefix}~{base}.somalier
    >>>

    output {
        File fingerprint = "~{prefix}~{base}.somalier"
    }
}

task CheckFingerprint {
    input {
        File fingerprint1
        File fingerprint2
        
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
        disk_gb: ceil(size(fingerprint1, "GB")+size(fingerprint2, "GB")) + 20,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "brentp/somalier:v0.2.19"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    command <<<
        set -euxo pipefail

        mkdir extracted
        mv ~{fingerprint1} extracted/
        mv ~{fingerprint2} extracted/

        somalier relate extracted/*.somalier
        tail -n1 somalier.pairs.tsv|cut -f3 > relatedness.txt
    >>>

    output {
        File relatednessfile = "somalier.pairs.tsv"
        String relatedness = read_string("relatedness.txt")
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
