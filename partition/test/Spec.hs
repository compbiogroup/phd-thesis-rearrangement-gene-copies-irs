import Control.Monad (unless)
import GenomesCheck as GC
import PartitionCheck as PC
import PFptCheck as FC
import System.Exit (exitFailure)

main :: IO ()
main = do
  ans <- fmap and . sequence $
    [ return True
    , GC.tests
    , PC.tests
    , FC.tests
    ]
  print ans
  unless ans exitFailure
