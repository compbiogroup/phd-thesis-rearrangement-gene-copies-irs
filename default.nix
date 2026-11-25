# Shell with everything needed to build and run the code of this repository:
#   nix-shell            (in the root of the repository)
# - cyc/: cmake, a C++ compiler with OpenMP (gcc, as in cyc/shell.nix) and make;
# - partition/: stack (the stack.yaml uses the Nix integration to get the GHC of
#   its resolver, so no GHC is needed here);
# - scripts/ and cyc/external/: bash, Python 3 with pandas (table_part_tese.py)
#   and coreutils (timeout, nproc);
# - db/: git-lfs, to download the data set.
{ pkgs ? import <nixpkgs> { } }:

pkgs.mkShell {
  name = "code-tese";
  nativeBuildInputs = with pkgs; [
    cmake
    gnumake
    gcc
    stack
    (python3.withPackages (ps: [ ps.pandas ]))
    bash
    coreutils
    git
    git-lfs
  ];
  shellHook = ''
    export CMAKE_CXX_COMPILER=${pkgs.gcc}/bin/c++
  '';
}
