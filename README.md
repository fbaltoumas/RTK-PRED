## RTK-PRED: automated detection and annotation of Receptor Tyrosine Kinases (RTKs) with profile Hidden Markov Models (pHMMs)

Online version: http://bioinformatics.biol.uoa.gr/RTK-PRED/

If you use RTK-PRED in your research please cite the following publication:

Filis, G., Baltoumas, F.A., Spanogiannis, G., Litou, Z.I., Iconomidou, V.A. (2023) **Proteome-Wide Detection and Annotation of Receptor Tyrosine Kinases (RTKs): RTK-PRED and the TyReK Database**. *Biomolecules*, 2023, 13(2), 270; [https://doi.org/10.3390/biom13020270](https://doi.org/10.3390/biom13020270)

Developers:
- Dr. Fotis Baltoumas (baltoumas@fleming.gr)
- Georgios Filis, Msc (georgefil28@gmail.com)
- Georgios Spanogiannis, Bsc (georgespano17@gmail.com)


#### Requirements
* **Operating System:** Linux, Microsoft Windows (native, with Windows Subsystem for Linux (WSL) or Cygwin), Apple Mac OS. BSD and other Unix-based systems will probably work too, in a manner identical (or at least very similar) to Linux and Mac OS.
* **Python 3** version **3.5** or higher
* **HMMER** version **3.1b1** or higher
* (*optional*) **Phobius** (alternatively, you can use the web version of Phobius using `-wp` or `--webphobius`, see below)

----

### First-time setup

#### 1. Download RTK-PRED
Download a zipped archive of the repository using the GitHub download button, or clone it with git:

    git clone https://github.com/fbaltoumas/RTK-PRED.git

#### 2. Install HMMER 3.1 or newer
Download and install [HMMER](http://hmmer.org/documentation.html).

**Note 1**: Older HMMER versions up to (and including) 3.0 DO NOT WORK, because they used different command line options.

**Note 2**: Windows users that do not use WSL need to compile HMMER from source.  However, you can find a pre-compiled Windows version of HMMER 3.2 [here](https://github.com/fbaltoumas/tools-for-windows/). This can be used with the Windows command line, Powershell or with CygWin.


#### 3. (optional) Install Phobius
Obtain a license, download and install Phobius from [here](https://phobius.sbc.su.se/data.html).  This step is optional, as RTK-PRED can also work with the Phobius web server as well, by using the `-wp` or `--webphobius` options.  However, running RTK-PRED with a local installation of Phobius is expected to be faster, especially for large datasets.

**Note**: Phobius is not currently available for native Windows or Mac OS.  Windows users can either use Windows Subsystem for Linux (Phobius will work with that) or, alternatively, use RTK-PRED with the web server option. Mac OS users will have to use RTK-PRED with the web server option.

#### 4. Configure RTK-PRED
Open `config.py` with a text editor and replace the following two paths:

    hmmer_dir = "/usr/bin/"
    phobius_dir = "./phobius/"

with the directories of HMMER and/or PHOBIUS in your machine, e.g.

    hmmer_dir = "/home/user/hmmer/bin"
    phobius_dir = "/home/user/phobius/"

You can also set the url for the Phobius web server (for use with the `-wp` version):

    phobius_url = "https://phobius.sbc.su.se/cgi-bin/predict.pl"

This can be useful, should the original url of the server change for any reason.

Alternatively, use the `--hmmerdir`, `--phobiusdir`, `-wp` or `--webphobius` options, to override `config.py`.

**Note:** If you don't have a Phobius license, you can still use RTK-PRED with the web-server version of Phobius through the `-wp` or `--webphobius` option.

----

### How to run:
USE:

    rtk-pred.py [-h] -i INPUT -o OUTPUT [--mkdir] [-wp] [--hmmerdir HMMERDIR] [--phobiusdir PHOBIUSDIR] [-v]

**required arguments:**

      -i INPUT, --input INPUT		input file in FASTA format
      -o OUTPUT, --output OUTPUT	the RTK-PRED output directory prefix

**optional arguments:**

      -h, --help            	show this help message and exit
      --mkdir               	if set, it automatically creates the output directory specified by '-o'.
      -wp, --webphobius     	use the web-server version of Phobius instead of a
	                        locally installed one (useful for Windows and MacOS
	                        machines, for which the 'decodeanhmm' binary used by
	                        Phobius is not available, as well as those users
	                        without a Phobius license), requires internet
	                        connection
      --hmmerdir HMMERDIR   	set the location of the HMMER binaries directory
                                (overrides config.py)
      --phobiusdir PHOBIUSDIR	set the location of the Phobius binaries directory
                                (overrides config.py), not applicable when
                                -wp/--webphobius is set. Also, not applicable in
                                Windows.
      -v, --version             display version and exit

The output files will consist of the following:

 1. summary.txt: The final RTK-PRED output, containing a summary of all predictions. **This is the file you most likely want to open and parse**.
 2. PTKhmm.res: HMMER results for the PTK pHMM, as run on the input.
 3. PTKs.fasta: FASTA file with all identified PTKs (single-domain PTKs, multi-domain PTKs and RTKs)
 4. singlePTKs.fasta: a subset of PTKs.fasta. FASTA file with PTKs containing only one PTK domain (single-domain PTKs and RTKs)
 5. phobius.txt: results from the Phobius transmembrane predictor, as run on singlePTKs.fasta
 6. RTKs.fasta: FASTA file containing all predicted RTKs.
 7. EC.res and JM.res: HMMER results for the extracellular domain (EC) and juxtamembrane region pHMMs (JM), as run on RTKs.fasta

 The final output (`summary.txt`) of RTK-PRED is in the following format:

    >> Query: sp|P00533|EGFR_HUMAN
    Classification: Receptor Tyrosine Kinase (RTK)
    SIGNALP:				From:1	To:24
    EC Domains (PF01030):	From:57	To:167	Recep_L_domain		Score: 105.2
    EC Domains (PF00757):	From:185	To:338	Furin-like		Score: 101.3
    EC Domains (PF01030):	From:361	To:480	Recep_L_domain		Score: 96.5
    EC Domains (PF00757):	From:474	To:557	Furin-like		Score: 23.8
    EC Domains (PF00757):	From:554	To:598	Furin-like		Score: 9.5
    TRANSMEM:				From:646	To:667
    Kinase Domain:				From:711	To:981		Score: 301.3
    RTK subfamily: Type 1 (EGF receptor subfamily)
    //
    >> Query: sp|P12931|SRC_HUMAN
    Classification: Non-Receptor Tyrosine Kinase (nRTK)
    Kinase Domain:				From:258	To:527		Score: 332.5
    //
    >> Query: sp|P04637|P53_HUMAN
    Classification: Not Tyrosine Kinase
    //
In each line you can see:

  1. Sequence Name
  2. Classification between RTKs, nRTKs, Not Tyrosine Kinase
  3. Tyrosine Kinase Domain Prediction (for RTKs and all PTKs)
  4. Tranmembrane Region Prediction (for RTKs only)
  5. Extracellular Domain Predictions (for RTKs only - provided there are positive hits)
  6. Classification between the 18 mammalian RTKs' subfamilies (for RTKs only - provided there are positive hits)

Each result starts with double greater-than sign `>>` and ends with double backslash sign `//`.

----

### Example Runs:

#### Use in Linux (all dependencies met):
Perform a run using the default installation of HMMER binaries (`/usr/bin/`) and a local installation of Phobius (e.g. `/home/user/phobius/`), defined in `config.py` as shown above:

    chmod +x rtk-pred.py
    ./rtk-pred.py -i test.fasta -o output --mkdir

#### Use in Linux (set custom paths):
Perform a run using by manually defining the HMMER and Phobius paths:

    chmod +x rtk-pred.py
    ./rtk-pred.py -i test.fasta -o output --mkdir --hmmerdir /home/user/hmmer/bin/ --phobiusdir /home/user/phobius/


#### Use in Linux with the Phobius web server:
Perform a run calling the web-server edition of Phobius (use of the `-wp` or `--webphobius` option):

    chmod +x rtk-pred.py
    ./rtk-pred.py -i test.fasta -o output --mkdir --hmmerdir /home/user/hmmer/bin/ -wp

#### Use in native Windows:
Open either the command line (CMD) or PowerShell and type the following:

    python3.exe rtk-pred.py -i test.fasta -o output --hmmerdir  "C:\Users\User\hmmer_windows\bin\" -wp --mkdir

where in `--hmmerdir` you define the location of your HMMER compiled files. Alternatively, you can edit `config.py` , enter the HMMER location there, and skip `--hmmerdir`.

**Note 1:** To use RTK-PRED in native Windows, you need HMMER compiled for Windows. You can find a Windows version of HMMER 3.2 in [this repository](https://github.com/fbaltoumas/tools-for-windows/), or you can download the source code from hmmer.org and compile it yourself.

**Note 2:** Since Phobius is only offered as a pre-compiled package for the Linux kernel, there is no Windows version available.  Therefore, you need to use the Phobius web server, with the option `-wp` or `--webphobius`.

**Note 3:** It is assumed that Python exists in your system path.  To add Python to the system path, either make sure you tick the `Add Python to Path` option during the installation of Python 3, or follow the procedure to add the Python3.exe executable to the system path (e.g., follow [this guide](https://answers.microsoft.com/en-us/windows/forum/windows_10-other_settings-winpc/adding-path-variable/97300613-20cb-4d85-8d0e-cc9d3549ba23)).

**Note 4:** In `--hmmerdir` (or in `config.py`), make sure you use Windows-style paths ("\\" instead of "/").

#### Use in Windows with Cygwin:
Open Cygwin and type:

       chmod +x rtk-pred.py
       ./rtk-pred.py -i test.fasta -o output --hmmerdir  /cygdrive/c/Users/User/hmmer_windows/bin/ -wp --mkdir

where in `--hmmerdir` you define the location of your HMMER compiled files. Alternatively, you can edit `config.py` , enter the HMMER location there, and skip `--hmmerdir`.

**Note 1:** To use RTK-PRED in Cygwin, you need HMMER compiled for Windows. You can find a Windows version of HMMER 3.2 in [this repository](https://github.com/fbaltoumas/tools-for-windows/), or you can download the source code from hmmer.org and compile it yourself.

**Note 2:** There is no Cygwin version available for Phobius.  Therefore, you need to use the Phobius web server, with the option `-wp` or `--webphobius`.

**Note 3:** By default, Cygwin does not install Python. Make sure that you have Python 3 installed.

**Note 4:** In `--hmmerdir` (or in `config.py`), make sure you use Unix-style paths (`/cygdrive/c/...` etc for Widnows paths).

#### Use in Windows with Windows Subsystem for Linux (WSL):
The WSL overlay essentially builds a Linux compatibility layer, with the native Linux kernel.  This means that binaries compiled for Linux run in WSL natively.  Therefore, all of the runs described in the **Use in Linux** subsections above also work with WSL in the same manner as an actual Linux OS. Most importantly, the local version of Phobius runs on WSL installations with no problems; therefore, you can use RTK-PRED without needing the web server.

#### Use in Mac OS:
Perform a run calling the web-server edition of Phobius (use of the `-wp` or `--webphobius` option):

    chmod +x rtk-pred.py
    ./rtk-pred.py -i test.fasta -o output --mkdir --hmmerdir /home/user/hmmer/bin/ -wp

**Note:** There is no Mac OS version available for Phobius.  Therefore, you need to use the Phobius web server, with the option `-wp` or `--webphobius`.
