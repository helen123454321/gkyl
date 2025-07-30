--------------------------------------------------------------------------------
-- App dependencies
--------------------------------------------------------------------------------
-- Load the App and equations to be used
local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler
-- local const = require "Lib.Constants"

local Te_Ti = 1.0 -- ratio of electron to ion temperaute
--local machNum = 1.5 -- Mach number computed from ion thermal speed
local n0 = 1.0 -- initial number density

local vd_e, vd_i = 0.0, 0.0

-- physical parameters (commented out are the SI constants if need be)
local me, qe = 1.0, -1.0 -- const.ELECTRON_MASS, -const.ELEMENTARY_CHARGE
local Te = 1.0 -- 1.5*const.ELEMENTARY_CHARGE
local vthe = math.sqrt(Te/me)

local mi, qi = 1836.0, 1.0 -- 39.948*const.MASS_UNIT, const.ELEMENTARY_CHARGE -- argon
local Ti = 1.0
local vthi = math.sqrt(Ti/mi)

local epsilon0 = 1.0 -- const.EPSILON0
local mu0 = 1.0
local wpe = math.sqrt(qe^2*n0/(epsilon0*me))
local lambdaD = vthe/wpe
local L = 128.0*lambdaD*2
local Nx = 128*2
local tEnd = 2000.0/wpe
local nFrame = 10

local gasGamma = 5.0/3.0
-- A little hack for the BCs
local function getAbsorbFunc(n, e)
   return function (dir, tm, idxIn, fIn, fOut)
      fOut[1] = n
      fOut[2] = 0.0
      fOut[3] = 0.0
      fOut[4] = 0.0
      fOut[5] = e
   end
end

-- Load data for initial Robertson profiles
-- I removed this 7/30/25 <3

--------------------------------------------------------------------------------
-- App construction
--------------------------------------------------------------------------------
local app = Moments.App.new {
   logToFile = false,
   tEnd = tEnd,
   nFrame = nFrame,
   lower = {0},
   upper = {L},
   cells = {Nx},
   timeStepper = "fvDimSplit",
   periodicDirs = {},
   
   elc = Moments.Species.new {
      charge = qe, mass = me,
      equation = Euler.new { gasGamma = gasGamma },
      equationInv = Euler.new { gasGamma = gasGamma, numericalFlux = "lax" },
      init = function (t, z)
         local rhoe = n0 * me
         local e = n0 * Te / (gasGamma - 1.0)
         return rhoe, 0.0, 0.0, 0.0, e
      end,

      evolve = true, 
      bcx = { Moments.Species.bcCopy, { getAbsorbFunc(n0*me*1e-10, n0*Te/(gasGamma-1.0)*1e-10) } },
   },

   ion = Moments.Species.new {
      charge = qi, mass = mi,
      equation = Euler.new { gasGamma = gasGamma },
      equationInv = Euler.new { gasGamma = gasGamma, numericalFlux = "lax" },
      init = function (t, z)
         local rhoi = n0 * mi
         local rhoiui = rhoi * vd_i -- vd_i = 0 by default
         local e = n0 * Ti / (gasGamma - 1.0)
         return rhoi, rhoiui, 0.0, 0.0, e
      end,
      evolve = true,
      bcx = { Moments.Species.bcCopy, { getAbsorbFunc(n0*mi*1e-10, n0*Ti/(gasGamma-1.0)*1e-10) } },
   },

   -----------------------------------------------------------------------------
   -- Electromagnetic field
   -----------------------------------------------------------------------------
   field = Moments.Field.new {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, z)
         return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      end,

      evolve = true,
      bcx = { Moments.Field.bcOpen, Moments.Field.bcReflect },
   },

}

app:run()