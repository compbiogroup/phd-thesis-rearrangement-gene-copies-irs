let
  pkgs = import <nixpkgs> { };
in
pkgs.mkShell {
  name = "haskell stack shell";
  buildInputs = with pkgs; [
    stack
    ghc
    haskell-language-server
    # haskell.compiler.ghc963
    # (haskell-language-server.override { supportedGhcVersions = [ "963" ]; })
  ];
}
