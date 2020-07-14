#!/usr/bin/env python
# -*- coding: utf-8 -*-

####### 1. Train Koopman operator Kd from simulation data

import KoopmanTraining

####### 2. Measure prediction error using trained Koopman operator Kd

import PredictionError

####### 3. Compare error bounds against actual error as a function of time

import Plot_ErrorBounds
