#ifndef UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H
#define UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H


#include "lib_disc/function_spaces/integrate.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "../../XBraidForUG4/src/interface/spatial_norm.h"
#include "../../XBraidForUG4/src/util/parallel_logger.h"

#include "biot_error_data.h"


namespace ug { namespace XBraidPoroelasticity {

    template<typename TDomain, typename TAlgebra>
    class BiotBraidSpatialNorm : public XBraidForUG4::BraidSpatialNorm<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        //--------------------------------------------------------------------------------------------------------------

        int m_uorder = 4;
        int m_porder = 2;
        double p_factor = 1;
        double u_factor = 1;


        //--------------------------------------------------------------------------------------------------------------

        BiotBraidSpatialNorm() : XBraidForUG4::BraidSpatialNorm<TDomain, TAlgebra>() {}

        ~BiotBraidSpatialNorm() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_order(int uorder, int porder) {
            this->m_uorder = uorder;
            this->m_porder = porder;
        }

        void set_parameter(double alpha, double lambda, double mu) {
            p_factor = alpha;
            u_factor = lambda + 2 * mu;
        }

        double norm(SP_GridFunction u) override {
            double norm_x = H1SemiNorm(*u.get(), "ux", this->m_uorder);
            double norm_y = H1SemiNorm(*u.get(), "uy", this->m_uorder);
            double norm_p = L2Norm(*u.get(), "p", this->m_porder);


            double pnorm = p_factor * norm_p * norm_p; //p_factor*unorm_p*unorm_p;
            double unorm = u_factor * (norm_x * norm_x + norm_y * norm_y);

            double total_norm = sqrt(unorm + pnorm);

            return total_norm;
        }
        //--------------------------------------------------------------------------------------------------------------
    };
}}

#endif //UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_ESTIMATOR_H
