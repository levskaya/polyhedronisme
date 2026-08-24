'use strict'
const pegjs = require('pegjs')

module.exports = function (source) {
  if (this.cacheable) this.cacheable()
  return pegjs.generate(source, {
    format: 'commonjs',
    output: 'source',
    optimize: 'speed'
  })
}
