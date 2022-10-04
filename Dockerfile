FROM continuumio/miniconda3:latest

WORKDIR /app

# Install create_aircraft requirements
COPY environment.yml /app/environment.yml
RUN conda config --add channels conda-forge \
    && conda env create -n create_aircraft -f environment.yml \
    && rm -rf /opt/conda/pkgs/*

# Copy all files after to avoid rebuild the conda env each time
COPY . /app/

# activate the create_aircraft environment
ENV PATH /opt/conda/envs/create_aircraft/bin:$PATH

# Launch the API

CMD ["bash", "-c", "source activate create_aircraft && flask run"]