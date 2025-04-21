module JuChrom

using Reexport

using Pkg.Artifacts
import BasicInterpolators
import Optim
import BSplineKit
import BSplineKit: diff
import Roots

@reexport using Unitful

const agilent = artifact"Agilent"
const andi = artifact"ANDI"
const calibration = artifact"calibration"
const shimadzu = artifact"Shimadzu"

end