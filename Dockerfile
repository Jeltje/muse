FROM    ubuntu

RUN     apt-get update && apt-get install -y wget zlib1g-dev

RUN     mkdir /opt/bin
RUN     wget -O /opt/bin/MuSEv1.0rc http://bioinformatics.mdanderson.org/Software/MuSE/MuSEv1.0rc_submission_c039ffa
RUN     chmod +x /opt/bin/MuSEv1.0rc

ENV     PATH /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/bin

# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* 

ENTRYPOINT ["/opt/bin/MuSEv1.0rc"]

