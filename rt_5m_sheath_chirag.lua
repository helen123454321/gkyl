--------------------------------------------------------------------------------
-- App dependencies
--------------------------------------------------------------------------------
-- Load the App and equations to be used
local Moments = G0.Moments
local Euler = G0,Moments
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
local fh = io.open("ni.txt")
local n_i, i = {}, 1
for l in fh:lines() do
   n_i[i] = l
   i = i+1
end

local fh = io.open("ne.txt")
local n_e, i = {}, 1
for l in fh:lines() do
   n_e[i] = l
   i = i+1
end

local fh = io.open("E.txt")
local E0, i = {}, 1
for l in fh:lines() do
   E0[i] = l
   i = i+1
end

local fh = io.open("ui.txt")
local ui, i = {}, 1
for l in fh:lines() do
   ui[i] = l
   i = i+1
end

--------------------------------------------------------------------------------
-- App construction
--------------------------------------------------------------------------------
local app = Moments.App {
   logToFile = false,
   tEnd = tEnd,
   nFrame = nFrame,
   lower = {0},
   upper = {L},
   cells = {Nx},
   timeStepper = "fvDimSplit",
   periodicDirs = {},
   
   elc = Moments.Species {
      charge = qe, mass = me,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      init = function (t, z)
         local x = z[1]
         local idx = x + 0.5 - 128
         if idx < 1 then
            idx = 1
         end
         if idx > 128 then
            idx = 128
         end         
         local rhoe = n_e[idx] * me         
         local e = n0*Te/(gasGamma-1.0)
         return rhoe, 0.0, 0.0, 0.0, e
      end,
      evolve = true, 
      bcx = { Moments.Species.bcCopy, { getAbsorbFunc(n0*me*1e-10, n0*Te/(gasGamma-1.0)*1e-10) } },
   },

   ion = Moments.Species {
      charge = qi, mass = mi,
      equation = Euler { gasGamma = gasGamma },
      equationInv = Euler { gasGamma = gasGamma, numericalFlux = "lax" },
      init = function (t, z)
         local x = z[1]
         local idx = x + 0.5 - 128
         if idx < 1 then
            idx = 1
         end
         if idx > 128 then
            idx = 128
         end
         local rhoi = n_i[idx] * mi         
         local rhoiui = rhoi * ui[idx]          
         local e = n0*Ti/(gasGamma-1.0)
         return rhoi, rhoiui, 0.0, 0.0, e
      end,
      evolve = true,
      bcx = { Moments.Species.bcCopy, { getAbsorbFunc(n0*mi*1e-10, n0*Ti/(gasGamma-1.0)*1e-10) } },
   },

   -----------------------------------------------------------------------------
   -- Electromagnetic field
   -----------------------------------------------------------------------------
   field = Moments.Field {
      epsilon0 = epsilon0, mu0 = mu0,
      init = function (t, z)
        local x = z[1]
        local idx = x + 0.5 - 128
        if idx < 1 then
           idx = 1
        end    
        local Ex = E0[idx] * 1.0        
         return Ex, 0.0, 0.0, 0.0, 0.0, 0.0
      end,
      evolve = true,
      bcx = { Moments.Field.bcOpen, Moments.Field.bcReflect },
   },

   emSource = Moments.CollisionlessEmSource {
      species = {"elc", "ion"},
      timeStepper = "direct",
   },
}

app:run()
