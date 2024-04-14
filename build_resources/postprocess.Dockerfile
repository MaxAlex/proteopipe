FROM python:3.9

RUN apt update && apt install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

RUN python -m pip install pandas scikit-learn numpy mokapot==0.9.1 intervaltree pyteomics

COPY run_mokapot.py /
COPY mod_mass.py /
COPY dinosaur_quantitation.py /
COPY comet.params.new /
COPY comet_search.py /
COPY call_pipeline_func.py /
COPY generate_final_report.py /
COPY isobaricquant_format.py /
COPY peptide_protein_aggregation.py /

COPY extract_isobaric_quant.py /
COPY integrate_isobaric_data.py /
