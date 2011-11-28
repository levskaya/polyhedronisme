
limiter = 999
t = undefined
verts = []
P = undefined

#vertices = []
diagonals = []


XOR = (a, b) -> (a or b) and not (a and b)

Area2 = (a,b,c) -> (b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])
Left      = (a, b, c) ->  Area2(a, b, c) > 0
LeftOn    = (a, b, c) ->  Area2(a, b, c) >= 0
Collinear = (a, b, c) ->  Area2(a, b, c) is 0

InCone = (a, b) -> #idxfunc
  a1 = (a+1+facelen)%facelen
  a0 = (a-1+facelen)%facelen
  if LeftOn(verts[a], verts[a1], verts[a0])
    return (Left(verts[a], verts[b], verts[a0]) and Left(verts[b], verts[a], verts[a1]))
  not (LeftOn(verts[a], verts[b], verts[a1]) and LeftOn(verts[b], verts[a], verts[a0]))

Between = (a, b, c) -> #vecfunc
  return false if Collinear(a, b, c)
  unless a[0] is b[0]
    (a[0] <= c[0]) and (c[0] <= b[0]) or (a[0] >= c[0]) and (c[0] >= b[0])
  else
    (a[1] <= c[1]) and (c[1] <= b[1]) or (a[1] >= c[1]) and (c[1] >= b[1])

IntersectProp = (a, b, c, d) -> #vecfunc
  return false if Collinear(a, b, c) or Collinear(a, b, d) or Collinear(c, d, a) or Collinear(c, d, b)
  XOR(Left(a, b, c), Left(a, b, d)) and XOR(Left(c, d, a), Left(c, d, b))

Intersect = (a, b, c, d) ->
  if IntersectProp(a, b, c, d)
    true
  else
    if Between(a, b, c) or Between(a, b, d) or Between(c, d, a) or Between(c, d, b)
      true
    else
      false

Diagonalie = (a, b) -> #idxfunc
  c = 0
  loop
    c1 = (c+1+facelen)%facelen
    if (c isnt a) and (c1 isnt a) and (c isnt b) and (c1 isnt b) and IntersectProp(verts[a], verts[b], verts[c], verts[c1])
      return false
    c  = (c+1+facelen)%facelen
    break unless c isnt 0
  true

Diagonal = (a, b) -> InCone(a, b) and InCone(b, a) and Diagonalie(a, b)

EarInit = ->
  v1 = vertices[0]
  loop
    v2 = v1.next
    v0 = v1.prev
    v1.ear = Diagonal(v0, v2)
    v1 = v1.next
    break unless v1 isnt verts[0]

Triangulate = ->
  n = vertices.length
  EarInit()
  z = limiter
  while z > 0 and n > 3
    z -= 1
    v2 = head
    y = limiter
    loop
      y -= 1
      broke = false
      if v2.ear
        v3 = v2.next
        v4 = v3.next
        v1 = v2.prev
        v0 = v1.prev
        diagonals.push [ v1, v3 ]
        v1.ear = Diagonal(v0, v3)
        v3.ear = Diagonal(v1, v4)
        v1.next = v3
        v3.prev = v1
        head = v3
        n--
        broke = true
      v2 = v2.next
      break unless y > 0 and not broke and v2 isnt head






