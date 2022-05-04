FROM ubuntu:18.04

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get install -y curl \
    && rm -rf /var/lib/apt/lists/*

COPY --from=quay.io/iobio/htslib:1.11 /bcftools-1.11/bcftools /usr/bin/bcftools

RUN curl -L https://github.com/brentp/vcfanno/releases/download/v0.3.3/vcfanno_linux64 -o /usr/bin/vcfanno \
    && chmod +x /usr/bin/vcfanno

RUN curl -L https://github.com/brentp/slivar/releases/download/v0.2.7/slivar -o /usr/bin/slivar_static \
    && chmod +x /usr/bin/slivar_static

COPY annotate_dev.sh /
COPY annotate_prod.sh /
