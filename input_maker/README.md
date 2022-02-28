# An input maker for Gaussian and ORCA calculations

This tool can be used to automate the extraction of coordinates from the outputs of two major quantum chemical software packages and create input files for the other steps of the workflow. It can also manipulate input files to change command line options or to convert between the formats of Gaussian and ORCA. Note that some functions are still work in progress.

## Requirements

The tool requires a basic Python installation (tested on 3.7+).

## Usage

- One functionality is to extract the last geometry from a Gaussian output file and write a basic Gaussian input (.gjf), ORCA (.inp) input or .xyz file. The created file will inherit the name of the Gaussian output file, appended with "_coords". For example, a file called "output_coords.inp" will be created if the first line of the following is entered:

```
makeinput -i output.log -o    # to obtain ORCA input
makeinput -i output.log -g    # to obtain Gaussian input
makeinput -i output.log -xyz  # to obtain an xyz file
```

- If no output file is pecified through the -i option, all the files with .log extension in the folder will be processed.

```
makeinput -o
makeinput -g
makeinput -xyz
```

- Another functionality takes a Gaussian input/output file together with a sample file that contains the command line and inserts the extracted charge, multiplicity and coordinates into the sample. There are also placeholder strings for solvent, gen basis (for Gaussian inputs) and external .xyz cordinates (for ORCA). Some examples can be found under the samples folder. This functionality can be used through the following:

```
makeinput -i input.gjf -s sample.gjf   # or -s sample.inp
makeinput -i input.xyz -s sample.gjf   # or -s sample.inp
makeinput -i output.log -s sample.gjf  # or -s sample.inp
```

> **_NOTE:_**  Solvents for the command line can be inserted using a placeholder 'SLVT'. The way it works is that a dictionary containing {'molecule name': 'solvent'} is defined manually in the *get_solvent()* function of *File_operat_mi.py*. If *'molecule name'* is found in the name of the file specified after -i, the corresponding solvent name is inserted into the placeholder string.

- The sample file can be applied to all Gaussian input (.gjf) files in the folder:

```
makeinput -s sample.gjf    # or -s sample.inp
```



