#ifndef UG_PLUGIN_XBRAIDBIOT_BIOT_ERROR_DATA_H
#define UG_PLUGIN_XBRAIDBIOT_BIOT_ERROR_DATA_H

#include <math.h>

#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug { namespace XBraidPoroelasticity {

    template<typename TDomain, typename TAlgebra>
    class BiotErrorData {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        double l2_norm_p = 0.0;
        double l2_norm_ux = 0.0;
        double l2_norm_uy = 0.0;

        double h1_norm_ux = 0.0;
        double h1_norm_uy = 0.0;

        int porder = 2;
        int uorder = 4;

        //--------------------------------------------------------------------------------------------------------------

        void compute(SP_GridFunction u) {

            this->l2_norm_p = L2Norm(*u.get(), "p", porder);;
            this->l2_norm_ux = L2Norm(*u.get(), "ux", uorder);
            this->l2_norm_uy = L2Norm(*u.get(), "uy", uorder);

            this->h1_norm_ux = H1SemiNorm(*u.get(), "ux", uorder);
            this->h1_norm_uy = H1SemiNorm(*u.get(), "uy", uorder);

        }

        void set_order(int porder = 2, int uorder = 4) {
            this->porder = porder;
            this->uorder = uorder;
        }
        //--------------------------------------------------------------------------------------------------------------
    };

}}

#endif //UG_PLUGIN_XBRAIDBIOT_BIOT_ERROR_DATA_H
