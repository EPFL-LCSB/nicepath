NICEpath - find pathways in a biochemical network
---------------------------------------------------
NICEpath is a pathway search tool for biochemical networks.
To install the dependencies, go to the nicepath directory and run:
```
$ cd /path/to/nicepath
$ make
```

To run NICEpath, go to the ./nicepath/ folder and run:
```
$ cd nicepath
$ python main.py INPUT-NAME
```

To search for pathways, you can 
1) use the available reaction networks (KEGG, available Retrobiosynthesis networks) --> set INPUT-NAME to KEGG or Retrobio 
	
2) add your own biochemical network, including weighted reactant-product pairs, and set the parameters according to your search criteria in ./input/INPUT-NAME

**Input:**
1) Reaction network, in the format of a BNICE.ch systemFile.txt (in ROOT/PROJECT/COMPOUND/)
2) Starting compounds for pathway search : source_compounds.txt
3) Parameter settings: parameters.txt

**Output:** in ./data/

Structure: ROOT/PROJECT/COMPOUND

ROOT: location of your project (default: ./data), 
PROJECT: name of your project (default: Retrobio), 
COMPOUND: name of your target compound (default: 3HP, test case)

