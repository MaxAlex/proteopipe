FROM adoptopenjdk/openjdk11:jdk-11.0.8_10-ubuntu

RUN apt update && apt install -y curl unzip
RUN apt install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install


RUN curl -L "https://github.com/fickludd/dinosaur/releases/download/1.2.0/Dinosaur-1.2.0.free.jar" -o "Dinosaur-1.2.1.free.jar"
