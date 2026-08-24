'use strict'
const chalk = require('chalk')
const semver = require('semver')
const packageConfig = require('../package.json')

module.exports = function () {
  const warnings = []

  if (packageConfig.engines && packageConfig.engines.node) {
    const currentNode = semver.clean(process.version)
    if (!semver.satisfies(currentNode, packageConfig.engines.node)) {
      warnings.push('node: ' +
        chalk.red(currentNode) + ' should be ' +
        chalk.green(packageConfig.engines.node)
      )
    }
  }

  if (warnings.length) {
    console.log('')
    console.log(chalk.yellow('To use this template, you must update following to modules:'))
    console.log()

    for (let i = 0; i < warnings.length; i++) {
      console.log('  ' + warnings[i])
    }

    console.log()
    process.exit(1)
  }
}
