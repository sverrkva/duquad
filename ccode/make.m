% Makefile

mex -O CC=gcc LD=gcc COPTIMFLAGS='-O3 -DNDEBUG'...
    CFLAGS="\$CFLAGS -std=c99  -I src/include"...
    src/main.c...
    src/print.c...
    src/fgm.c...
    src/gdm.c...
    src/math_functions.c... 
    src/dgm.c...
    src/general_functions.c...
    src/dfgm.c...
    src/alm.c...
    src/falm.c
