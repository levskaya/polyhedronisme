# Polyhédronisme
#===================================================================================================
#
# A toy for constructing and manipulating polyhedra and other meshes
#
# Includes implementation of the conway polyhedral operators derived
# from code by mathematician and mathematical sculptor
# George W. Hart http://www.georgehart.com/
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License


# Testing Functions
#===================================================================================================

#report on face topology
topolog = (poly) ->
  str=""
  for f in poly.face
    for v in f
      str+="#{v}->"
    str+="\n"
  console.log str

testrig = ->
  seeds=["T","O","C","I","D","P3","P4","A4","A5","Y3","Y4"]
  ops = ["k","a","g","p","d","n","x","*"]
  console.log "===== Test Basic Ops ====="
  for op in ops
    console.log "Operator #{op}"
    for seed in seeds
      console.log op+seed+":", generatePoly op+seed
  console.log "===== Done Testing Basic Ops ====="
