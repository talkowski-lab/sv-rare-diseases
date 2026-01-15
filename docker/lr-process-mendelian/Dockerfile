FROM quay.io/ymostovoy/lr-utils-basic:2.0

COPY ./dup_to_ins_LR.py /
COPY ./ins_to_dup.py /
COPY ./ins_to_dup_collapse.py /
COPY ./symbalts.py /
COPY ./abs_svlen.py /
COPY ./orig_dup_to_dup.py /

RUN apt-get update && apt-get install python3-pip -y

RUN pip3 install pysam
