# Use an official R base image
FROM r-base:4.4.1

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git

# Install remotes package, SAIGE_META, and SKAT from GitHub
RUN R -e "install.packages('remotes', repos='http://cran.rstudio.com/'); \
           library(remotes); \
           install_github('leelabsg/SAIGE_META'); \
           install_github('leelabsg/SKAT')"

# Create directories for the files
RUN mkdir -p /app/R /app/inst/scripts

# Copy the specific files
COPY R/MetaSAIGE.R /app/R/
COPY inst/scripts/RV_meta_GC.R /app/inst/scripts/

# Set the default command
CMD ["bash"]
