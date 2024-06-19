[![DOI](https://zenodo.org/badge/678961833.svg)](https://zenodo.org/doi/10.5281/zenodo.12163217)
# Innovation Processes for Inference

This repository contains the code associated to [Innovation Processes for Inference (2023)](https://arxiv.org/abs/2306.05186) by G. Tani Raffaelli, M. Lalli, and F. Tria

Code from Giulio Tani Raffaelli

## Citation

If you use this code or the data we prepared for academic work please cite: *G. Tani Raffaelli, M. Lalli, and F. Tria, Innovation processes for inference (2023), arXiv:2306.05186*

## Installation

1. Clone the repository in a new directory
1. Navigate to the `InnovationProcessesInference` folder
1. Run `make`, this will compile the C/C++ code and create the two default directories `data` and `res` in the parent directory
    - **NOTE:** While data can reside anywhere, currently there is no way to specify a different path for the results files other than the default. In any case, for most usage cases you won't need to access that folder.
    - **NOTE:** to change the normalisation of $P_0$ open `InnovationProcessesInference/bookprob.cpp` and edit the line `#define _P0_NORMALIZATION_ 0`, then recompile.

### Package requirements

Last tested on Ubuntu 22.04
- Python>=3.9.16
- Numpy>=1.20
- Pandas>=1.5.2
- C++17 (recommended g++=11.4.0)

## Testing

1. Extract the archive `InnovationProcessesInference/sample_data/English_literature/englis_literature.zip` in the folder `data`
1. Assuming that you are in the same directory containing `InnovationProcessesInference`, `data` and `res` run:
    ``` bash
    export PYTHONPATH="./InnovationProcessesInference:$PYTHONPATH"
    python3 -m cp2d config
    ```
    Follow the instructions to build the configuration file agreeing to a global configuration file. The program tries to be considerate of other users of shared machines when there's no workload manager.
1. Run:
    ```bash
    python3 -m cp2d attribute data/english_literature -n 9 -s test -l 1 -V
    ```
    In this case we are using OSF N-grams of length 9 (`-n 9`) in leave-one-out mode (`-l 1`) and a suffix `test` for the results folder (`-s test`), for further information about the options run `python3 -m cp2d attribute -h`
1. You'll see some information reminding you of the value of some compile-time variables. E.g., form `base_experiment.cpp` some possible flags are `_USE_JOINT_PARAMS_` and `_FLEXIBLE_AUTHOR_SIZE_` that refer to an alternative version where also the reference authors are split in fragments, `_USE_BOOK_COUNTS_` to count the tokens only once for every document in which they appear, `_FULL_AUTHOR_PARAMETERS_` allows to use the parameters computed on all the known texts from the author even when attributing its own books.
    - **NOTE:** author number 0 is special, this number is reserved for texts of unknown author. If you want to optimise the parameters with n-fold cross validation and the assign the author 0 *using the same parameters for the author's processes*, it's the time to un-comment `#define _FULL_AUTHOR_PARAMETERS_` and recompile.
1. Then you'll see some progress indicators, information about the execution (as the number of text-author comparisons) and other technical more useful with larger datasets.
1. After completion (about one minute) you'll see a report giving different metrics about attribution. In our case, since we are using full documents each column will contain identical values.
    - *Try:* run with `-f 1000` to observe the effect of different ways to combine the probabilities of different fragments (you'll see two distinct values), add `-F 50000` to split also the author in fragments, now all the row will contain different values. This last run can take ~10 time longer.

### More testing and insights

1. The number "1.0" that appears in the report means that we are using the plain CP2D, i.e. the $\delta-$ CP2D with $\delta = 1$. To change the value of $\delta$ to some different value, e.g. 0.7, just add `-D 0.7` when calling the program.
    - **NOTE:** the intermediate results of the computation are cached. This will speed up considerably the execution whe only changing $\delta$, but can take a lot of disk space (*hundreds of GB for the blog corpus in the paper!*), you are in charge of cleaning unwanted files from the `res` subdirectories.
    - **NOTE:** there is an attempt to prevent messing with the cache, however changing compile flags, editing the corpus, ignoring warnings and being creative can easily lead to inconsistent results.
    - *Tip:* you can compute the results with different values of $\delta$ in a single run, it's enough to ad a space separated list of values like `-D 0.7 0.9 1.1 1.3`. This is more efficient both for humans and the machine.
1. Have a look in the `res` folder you'll find a folder named `OSFG-09N-L1O-test`, it means that was produced with OSF N-grams (`OSFG`) of length 9 (`09N`) in a leave-one-out (`L1O`) procedure, then there is the suffix we specified before.
    - *Tip:* use the suffix to register what (version of) dataset you're using, it's a good way to keep track of what's where and avoids relying entirely on the program ability to recognise errors.
1. Look into `OSFG-09N-L1O-test`, you'll see at least four things:
    1. `parameters.json`, for every run that ends up in the same folder adds a line summarising the parameters of the computation. It's used by the program to attempt to prevent errors but can be useful for humans too.
    1. `results.json`, for every run that ends up in the same folder adds a line summarising the results along with some parameters like the value of $\delta$, documents excluded or not assigned and more. If the combination of parameters when you launch the program matches a line in this file, the results are drawn from here, if in the meantime you changed compile flags or the corpus and ignore the warnings, the results will be wrong.
        - **NOTE:** instead of $\delta$ what's reported is $\log_{10}\delta$
        - **NOTE:** the `"unknowns"` element will contain the attributions of the documents with author number 0.
    1. `seq`, a folder containing a file for every author with all its texts partly preprocessed and in a format convenient for computing attributions (if the corpus is large the `seq` folders in the various folders of `res` can take up quite some space, clean what you don't need), and a binary file `slices.bin` with information for the program about which document is in which slice (in leave-one-out every document has its own slice), copy this around and use the `-M` flag to ensure that the folds in which the corpus is divided are the same across tries.
    1. `fullA`, a folder with all the results when the author is not split in fragments, there might be also another folder `050000BA` if you tried `-F 50000`. Lets look inside `fullA`, we'll see:
        1. `params.bin`, a binary file with the $\alpha$ and $\theta$ parameters for all authors in all slices,
        1. `infBF`, a folder containing the results when using full documents. There might be also another folder `001000BF` if you tried `-f 1000`. Let's look inside `infBF`, we'll see:
            1. `MR_Attribution.txt`, this file was produced by the flag `-V` and contains detailed information about which document was assignet to which author, how far is the second, ... If the corpus is large can be big and not very handy but in small corpora is useful to spot if some author is attracting all attributions.
            1. `OSFG-09N-L1O-test2-fullA-infBF_fra.res`, machine readable information about the number of fragments of every document.
            1. `OSFG-09N-L1O-test2-fullA-infBF_ids.res` and `OSFG-09N-L1O-test2-fullA-infBF_pro.res`, cached values of probabilities in an optimised format (do not expect to manage to compress the much more than 20% in a very slow way), used to recompute the results quickly if only $\delta$ changes.
            1. `params.bin`, same as above.
            1. `results.json`, a subset of the one seen above just for this combination of author and document fragments.
            1. `slices.bin`, the last slicing used for this combination of author and document fragments.
1. You may have noticed in the `"exc": [], "noa": [], "asso": {}` elements in the results file, `"asso"` will change if providing an "author association file" with the option `-a`: you can treat some authors as different during probability computation but then as the same author during attribution, this option is here for this. The other two are not accessible via command line but require the use of a python script.
    - *Tip:* You are curious to see if the episode of the death of Dido gets assigned to the correct book of the Aeneid or *at least* to the Aeneid at all? An "author association file" can come in handy. Use the name of a file contained in the same folder of the database (trying to helping in avoiding mess) which contains a json dictionary in the format {"old_auth_num":<new_auth_num>} where the old "authors" are books and the "new" will be the same for all the books of the Aeneid. (This is just an idea. There are other ways to answer the above question with different theoretical implications.)

### Programmatic approach

To deploy the full possibilities of CP2D import `cp2d` as a module, create a `cp2dExperiment` object, `run()` it and extract the `results()`.
