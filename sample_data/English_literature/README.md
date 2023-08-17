# English literature corpus

The current folder contains the English books used in [Innovation Processes for Inference (2023)](https://arxiv.org/abs/2306.05186) by G. Tani Raffaelli, M. Lalli, and F. Tria in a `.zip` archive.

## Filenames format

The format of the file names is: A{author number}B{book number}.txt

The file authorNames.dat contains the association between author numbers and author names, the file bookNames.dat contains the association between pairs (author number, book number) and book titles.

## Technical info

All the texts were encoded as Latin1 (ISO/IEC 8859-1:1998), occasional out-of-encoding characters translitterated, and then prepared using the `cp2d prepare` command for further processing.

_NB:_ all newlines are removed from the texts, when opening with a text editor be sure that it behaves with lines over 1 MB. As an example, currently Sublime Text and VS Code open these files smoothly.