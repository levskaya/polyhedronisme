# Polyhédronisme Cake Build Script
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License
#

fs     = require 'fs'
{exec} = require 'child_process'
#{spawn} = require 'child_process'

appFiles  = [
  'geo'            # math, geometry
  'polyhedron'     # core mesh functions
  'topo_operators' # flagset->polyhedron constructor and topological operators
  'geo_operators'  # geometrical refiners/operators
  'triangulate'    # 2d delauney triangulation routine
  'testing'        # random testing funcs
  'csg'            # csg union,intersect,subtract routines
  'canvas_ui'      # 2d renderer and simple UI
]
outFile = 'polyhedronisme.coffee'
jsFile = 'polyhedronisme.js'
jsMinFile = 'polyhedronisme.min.js'

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

task 'minify', 'Minify the resulting application file after build', ->
  exec "yuicompressor -o \"#{jsMinFile}\" #{jsFile}", (err, stdout, stderr) ->
    throw err if err
    console.log stdout + stderr

#option '-e', '--environment [ENVIRONMENT_NAME]', 'set the environment for `task:withDefaults`'
#task 'task:withDefaults', 'Description of task', (options) ->
#  options.environment or= 'production'

#task 'minify2', 'Minify the resulting application file after build', ->
#  yui = spawn 'yuicompressor', ['-o', 'polyhedronisme.min.js', 'polyhedronisme.js']
#  yui.stdout.on 'data', (data)-> console.log('stdout: ' + data)
#  yui.stderr.on 'data', (data)-> console.log('stderr: ' + data)
#  yui.on 'exit', (code)-> console.log('child process exited with code ' + code)