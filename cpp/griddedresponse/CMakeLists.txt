# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

install(FILES griddedresponse.h aartfaacgrid.h dishgrid.h mwagrid.h
              phasedarraygrid.h skamidgrid.h
        DESTINATION "include/${CMAKE_PROJECT_NAME}/griddedresponse")
install(
  FILES aartfaacgrid.h
        airygrid.h
        dishgrid.h
        griddedresponse.h
        mwagrid.h
        phasedarraygrid.h
        skamidgrid.h
  DESTINATION "include/${CMAKE_PROJECT_NAME}/griddedresponse")
