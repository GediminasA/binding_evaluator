rule provean_nr:
    output:
        "data/nr/2011-08/nr.pal"
    shell:
        """
        mkdir --parents $(dirname {output})
        (
            cd $(dirname {output})
            for NUMBER in $(seq 0 5)
            do
                wget ftp://ftp.jcvi.org/data/provean/nr_Aug_2011/nr.0$NUMBER.tar.gz
                wget ftp://ftp.jcvi.org/data/provean/nr_Aug_2011/nr.0$NUMBER.tar.gz.md5
                md5sum --check nr.0$NUMBER.tar.gz.md5
                tar -xf nr.0$NUMBER.tar.gz
            done
        )
        """
