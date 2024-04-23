FROM tensorflow/tensorflow:2.13.0-gpu

RUN apt install ca-certificates gnupg
RUN gpg --homedir /tmp --no-default-keyring --keyring /usr/share/keyrings/mono-official-archive-keyring.gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb [signed-by=/usr/share/keyrings/mono-official-archive-keyring.gpg] https://download.mono-project.com/repo/ubuntu stable-focal main" | tee /etc/apt/sources.list.d/mono-official-stable.list
RUN apt update
RUN apt install -y mono-devel=6.12\*


# Pythonnet: 3.0.1 (from PyPI)
# Note: pycparser must be installed before pythonnet can be built
RUN pip install pycparser \
  && pip install pythonnet==3.0.1

RUN apt install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

#==1.1.2
RUN python -m pip install pandas numpy boto3 deeplc
RUN python -m pip install pyteomics

RUN rm -rf /var/lib/apt/lists/* /tmp/*

RUN git clone https://github.com/BlaisProteomics/multiplierz.git --branch patches
RUN python -m pip install -e multiplierz

RUN python -m pip install peptdeep

# COPY peptdeep_generic_models.tar.gz /
# RUN tar -xvf peptdeep_generic_models.tar.gz
RUN python -c "from peptdeep.model.ms2 import pDeepModel"  # Causes models to be downloaded, if they weren't

RUN python -m pip install mokapot==0.9.1 intervaltree

COPY annotate_psms_via_models.py /
COPY mod_mass.py /
COPY run_mokapot.py /
COPY call_pipeline_func.py /
