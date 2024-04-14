FROM python:3.9

RUN apt update && apt install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

RUN python -m pip install pandas scikit-learn numpy

COPY search_output_reconciliation.py /
COPY call_pipeline_func.py /
