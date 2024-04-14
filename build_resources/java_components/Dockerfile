FROM adoptopenjdk/openjdk11:jdk-11.0.8_10-ubuntu

RUN apt update && apt install -y curl unzip
RUN apt install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

RUN apt install -y python3 python3-pip
RUN python3 -m pip install pandas

COPY Dinosaur-1.2.1.free.jar /
COPY IsobaricQuant_LINUX.jar /
