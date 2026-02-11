# Minimal TeX Live 2020 profile for Docker
selected_scheme scheme-basic
TEXDIR /usr/local/texlive/2020
TEXMFCONFIG ~/.texlive2020/texmf-config
TEXMFVAR ~/.texlive2020/texmf-var
binary_x86_64-linux 1

# Only install the minimal collections required to run tlmgr
collection-basic 1
collection-latex 1
collection-latexrecommended 1
