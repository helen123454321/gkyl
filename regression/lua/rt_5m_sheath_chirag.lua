--------------------------------------------------------------------------------
-- App dependencies
--------------------------------------------------------------------------------
-- Load the App and equations to be used
local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler
-- local const = require "Lib.Constants"

local Te_Ti = 10.0 -- ratio of electron to ion temperaute (up to date, paper, but he has it as 10. in og code)
--local machNum = 1.5 -- Mach number computed from ion thermal speed
local n0 = 1.0 -- initial number density (up to og petr)

local vd_e, vd_i = 0.0, 0.0 -- not used? 

-- physical parameters (commented out are the SI constants if need be)
local me, qe = 1.0, -1.0 -- const.ELECTRON_MASS (up to og petr), -const.ELEMENTARY_CHARGE (up to og petr)
local Te = 1.0 -- 1.5*const.ELEMENTARY_CHARGE (up to date)
local vthe = math.sqrt(Te/me) -- thermal spped of electron (up to og petr) - H

local mi, qi = 1836.1, 1.0 -- 39.948*const.MASS_UNIT, const.ELEMENTARY_CHARGE -- argon (up to og petr)
local Ti = Te/Te_Ti  -- Temp of Ion - H 
local vthi = math.sqrt(Ti/mi) -- Thermal speed of ions 

local epsilon0 = 1.0 -- const.EPSILON0 -- Permittivity of free space. (up to og petr)
local mu0 = 1.0/100 -- Permittivity of free space. (up to og petr)
local wpe = math.sqrt(qe^2*n0/(epsilon0*me)) -- -- Plasma frequency. (up to og petr)
local lambdaD = vthe/wpe  --(up to og petr)

-- domain size
local L = 128.0*lambdaD*2 -- domain length BUT ONLY ONE! -H
local Nx = 128*2 -- grid size but only in one dimention?  -H
local tEnd = 400/wpe -- Chirag used to be 2000.0/wpe, but chagecd to match paper -H
local nFrame = 40
local dx = L/Nx -- cell spaceing (method copied from og petr)

local gasGamma = 5.0/3.0
-- A little hack for the BCs (he means bondry conditions )
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
   cells = {dx},
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

   source =  {
      {
      kind = "collisionless-em",
      species = {"elc", "ion"},
      timeStepper = "direct",
      },
   },

}

app:run() 