FROM python:3.9

RUN apt update && apt install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

RUN apt install -y wget
RUN wget https://github.com/UWPR/Comet/releases/download/v2023.01.2/comet.linux.exe
RUN chmod +x comet.linux.exe

# COPY multiplierz-2.2.1-py3-none-any.whl /
# RUN python -m pip install multiplierz-2.2.1-py3-none-any.whl
RUN git clone https://github.com/BlaisProteomics/multiplierz.git --branch patches
RUN python -m pip install -e multiplierz

RUN python -m pip install pandas numpy boto3

COPY utilities/comet.params.new /
COPY comet_search.py /
