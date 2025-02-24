# leven.c v1.0.2  
Levenshtein distance calculator  

Athanassios Protopapas  
Institute for Language & Speech Processing, February 2011  
Uploaded to github in February 2025  

This is optimized C code that calculates mean Levenshtein distances for a set of target items relative to a reference lexicon.  
Item files must have been preprocessed to remove nonletter characters and convert to lowercase (for orthographic strings).  

Compile with: `gcc -O9 -o levencmd leven.c`

The code was based on the [algorithm from Wikipedia](http://en.wikipedia.org/wiki/Levenshtein_distance). Thanks to Tal Yarkoni for help and examples.  

Calculation of distances for large sets of items is very time consuming; this code is highly optimized and **a lot** faster than plain coding of the algorithm.  

The code was uploaded to github after many years to facilitate discoverability and accessibility. It has long been [available](http://speech.ilsp.gr/iplr/leven-c.rar) via the [ILSP](https://www.ilsp.gr/en/home-2/) [IPLR](http://speech.ilsp.gr/iplr/index.htm) [downloads page](http://speech.ilsp.gr/iplr/downloads.htm). An optimized Windows executable (rar-archived) can also be [downloaded](http://speech.ilsp.gr/iplr/leven.rar) from there.
