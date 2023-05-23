#--------------------------------------------------------------------------------------------------------------
# Abstract gravity theories

abstract type AbstractGravity end

abstract type GeneralRelativity <: AbstractGravity end
GR = GeneralRelativity

abstract type DamourEspositoFarese <: AbstractGravity end
DEF = DamourEspositoFarese

abstract type extDamourEspositoFarese <: AbstractGravity end
extDEF = extDamourEspositoFarese

abstract type MendesOrtiz <: AbstractGravity end
MO = MendesOrtiz

#--------------------------------------------------------------------------------------------------------------
# Physics for gravity theories

