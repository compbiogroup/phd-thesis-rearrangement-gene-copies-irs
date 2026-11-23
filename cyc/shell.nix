with import <nixpkgs> { };
stdenv.mkDerivation {
# clangStdenv.mkDerivation {
  name = "env";
  nativeBuildInputs = [ cmake ];
  buildInputs = [ llvmPackages.openmp ];
  shellHook = ''
    export CMAKE_CXX_COMPILER=${pkgs.gcc}/bin/c++
  '';
}
