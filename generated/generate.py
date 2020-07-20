#!/usr/bin/env python3
import sys
import os
import errno
import argparse
import importlib.util


cmdLineParser = argparse.ArgumentParser()
cmdLineParser.add_argument('--order')
cmdLineParser.add_argument('--arch')
cmdLineParser.add_argument('--yateto')
cmdLineParser.add_argument('--outputDir')
cmdLineParser.add_argument('--matrixDir')
cmdLineArgs = cmdLineParser.parse_args()

sys.path.append(cmdLineArgs.yateto)


from yateto import *
from yateto.ast.visitor import PrettyPrinter, FindTensors, PrintEquivalentSparsityPatterns
from yateto.codegen.code import Cpp
from predictor import PredictorGenerator, SourceGenerator

arch = useArchitectureIdentifiedBy(cmdLineArgs.arch)

g = Generator(arch)

predictor = PredictorGenerator(cmdLineArgs.order,cmdLineArgs.matrixDir)
predictor.add(g)

source = SourceGenerator(cmdLineArgs.order,cmdLineArgs.matrixDir)
source.add(g)

gemm_cfg = predictor.gemm_cfg(arch)
g.generate(outputDir=cmdLineArgs.outputDir, gemm_cfg=gemm_cfg)

for kernel in g.kernels():
    title = 'AST of {}'.format(kernel.name)
    print(title)
    print('=' * len(title))
    PrettyPrinter().visit(kernel.ast)
    print(' ')
