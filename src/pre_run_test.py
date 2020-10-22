#!/usr/bin/python3
import os

def main():
    print('Biomage ingest running.')

    if not os.environ.get('EXPERIMENT_NAME'):
        print('The EXPERIMENT_NAME environment variable must be set.')
        print("Try running as: EXPERIMENT_NAME=\"Sample experiment\" docker-compose up --build")
        exit(1)
    
    if os.listdir('/output'):
        print(os.listdir('/output'))
        print('The output directory is not empty. Please clear it out before proceeding.')
        exit(1)
    
    print('Folders seem okay, starting R pre-processing...')

if __name__ == '__main__':
    main()

#r$> map[which(!is.na(map))] 

"""
r$> map <- mapIds(org.Gg.eg.db, 
                         keys=rownames(t), 
                         column="SYMBOL", 
                         keytype="ENSEMBL", 
                         multiVals="first") 
"""