# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

install(
  FILES aartfaacpoint.h
        airypoint.h
        dishpoint.h
        mwapoint.h
        pointresponse.h
        phasedarraypoint.h
        skamidpoint.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}/pointresponse")
