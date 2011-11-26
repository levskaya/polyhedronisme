# Polyhédronisme Cake Build Script
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License
#

fs     = require 'fs'
{exec} = require 'child_process'

appFiles  = [
  'geo'
  'polyhedron'
  'topo_operators'
  'geo_operators'
  'testing'
  'canvas_ui'
]
outFile = 'polyhedronisme.coffee'

task 'build', 'Build single application file from source files', ->
  appContents = new Array remaining = appFiles.length
  for file, index in appFiles then do (file, index) ->
    fs.readFile "#{file}.coffee", 'utf8', (err, fileContents) ->
      throw err if err
      appContents[index] = fileContents
      process() if --remaining is 0

  process = ->
    fs.writeFile "#{outFile}", appContents.join('\n\n'), 'utf8', (err) ->
      throw err if err
      exec "coffee --compile #{outFile}", (err, stdout, stderr) ->
        throw err if err
        console.log stdout + stderr
        fs.unlink "#{outFile}", (err) ->
          throw err if err
          console.log 'Done.'

