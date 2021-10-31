
#include <pybind11/pybind11.h>

#include "spm_bsplins.h"
#include "spm_bsplinc.h"
#include "spm_conv_vol.h"
#include "spm_slice_vol.h"


PYBIND11_MODULE(_rtspm, m) {
    m.doc() = R"pbdoc(
        Python adaptation of SPM for OpenNFT project
        --------------------------------------------

        .. codeauthor:: SPM Toolbox Team, Nikita Davydov and OpenNFT Team

        .. currentmodule:: rtspm
        .. autosummary::
           :toctree: _generate

           bsplins
           bsplins_multi
           bsplinc
           conv_vol
           slice_vol
    )pbdoc";

    m.def("bsplins", &spm_bsplins);
    m.def("bsplins_multi", &spm_bsplins_multi);
    m.def("bsplinc", &spm_bsplinc);
    m.def("conv_vol", &spm_conv_vol);
    m.def("slice_vol", &spm_slice_vol);
}
