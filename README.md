# vocab_create

## What is vocab_create

vocab_create is a multi-threaded C++ program to collect the most-common substrings from a set of text files. It is designed to be reasonably fast and use relatively little memory.

When linking against a Boost.IOStreams build with decompression support, files can be automatically decompressed when searching for substrings. Aside from the presence of compressed files, which will be decompressed, the assumption is that all files in the provided directories are text files. Furthermore, the assumption is that all text is UTF-8 (although some effort to skip random other bytes is made).

This was inspired by reading about the [Byte Pair Encoding](https://en.wikipedia.org/wiki/Byte_pair_encoding) (BPE) used by [Natural language processing](https://en.wikipedia.org/wiki/Natural_language_processing) (NLP) techniques (including some of today's [large language models](https://en.wikipedia.org/wiki/Large_language_model) (LLMs)). However, this program does not compute a BPE. Instead, it uses a data structure known as a [count-min sketch](https://en.wikipedia.org/wiki/Count%E2%80%93min_sketch) (CMS) to directly determine the most-common substrings (up to some specified number of characters) along with their relative frequency.

## Usage

``bash
$ vocab_create --log-level debug vocab.tsv /opt/data/

$ vocab_create -h
Usage:
  -h [ --help ]                     Display this help message
  --filter-buffer-size arg (=32768) Filter buffer size for decompression
  --log-level arg                   Output log level
  --log-file arg                    Output log file
  --threads arg                     Number of threads
  -c [ --config ] arg               Configuration file
  --max-token-length arg (=10)      Maximum token length
  --max-tokens arg (=10000)         Maximum number of tokens to generate
  --cms-width-factor arg (=100)     Width of the count-min sketch, 
                                    multiplicative to the maximum number of 
                                    tokens
  --cms-depth arg (=3)              Depth of the count-min sketch
  --output-file arg                 Output file
  --input-directory arg             Input directory

The input directory is iterated recursively for input files. All regular files
are treated as text files, with automated decompression provided for files
with extensions: .bz2 (bzip2) .gz (gzip) .lzma (lzma) .zz (zlib)

Note:
Options that can be provided on the command line can also be provided in a
configuration file. The configuration file uses a line-based 'name = value' syntax.
``

The output file will look something like this:

``bash
$ cat vocab.tsv
_	1
e	0.9288
a	0.902373
i	0.781876
o	0.700804
n	0.6902
r	0.625616
t	0.565587
1	0.564749
s	0.555351
l	0.460784
0	0.442342
.	0.377126
3	0.352061
...
ge	0.0270703
atio	0.0269839
ation	0.0269018
ss	0.0266672
so	0.0263907
24	0.0263628
ai	0.0262487
ig	0.025936
ist	0.025761
ter	0.0256278
ki	0.025567
ho	0.0254289
bo	0.0250904
16	0.02482
iv	0.0248142
21	0.0248003
pe	0.0247993
ent	0.0247608
...
``

## Maturity and Stability

This program is in an "I just got something working" state. No interface or functional stability should be presumed at this point in time. There are no official releases.

## Support

This library is being developed as a hobby of the author. Bugs should be reported using the repository bug tracker. Contributions are welcome.

