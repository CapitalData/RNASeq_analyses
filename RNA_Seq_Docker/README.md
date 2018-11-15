---
title: "RNA_Seq_Docker_Readme"
author: "Kiersten Henderson"
date: "11/13/2018"
output: pdf_document
---


## Docker Instructions

The docker directory is called "Actual_RNA_Seq_Docker".

Before creating the docker, the genome reference file you will use and the raw data must be placed into the appropriate directories under the data directory in the docker.

Note that the raw data should not be unzipped.

The RunRSEM.sh script has been modified to work in this docker and is called RunRSEM_for_Docker.sh.

When I (KH, Oct 2018) added it, I got it to the point of executing so that it doesn't error. But I never had it finish RSEM -  the script RunRSEM_for_Docker.sh must most likely be adapted to get the analysis to run to completion


Notes: 

It takes more than 20 minutes to create this docker and it is big.
The RunRSEM_for_Docker.sh script does not execute automatically, you must do that after the docker builds.

To build/delete docker: 

1) cd into docker directory (the one that contains the dockerfile ie. Acutal_RNA_Seq_Docker)

2) Run this command: docker build -t <name_your_docker> .

3) Confirm docker build and learn your docker's assigned image id with: docker image ls

4) To remove image: docker rmi <image id>

5) Before removing image, you may first need to delete container: docker rm <container id>

6) To work within the docker CLI and run the RSEM script:

    docker run -it <name_your_docker> bash
    --from docker cli you can navigate through docker 
    --execute RunRSEM_for_Docker.sh with: bash RunRSEM_for_Docker.sh
    
    
7) Do note that every time you are going to run this analysis, you'll have to delete your docker and build a new one where you add in the 
new data files you want to work with. This means you'll want to move any analysis files you generate out of the docker directory for storage.










