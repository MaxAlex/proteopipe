FROM chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

RUN apt update
RUN apt install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install