#ifndef UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_CONTROL_H
#define UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_CONTROL_H

#include "braid_biot_precomputed.h"
#include "../../Poroelasticity/src/biot_tools.h"
#include "../../Poroelasticity/src/barry_mercer.h"

#include "../../XBraidForUG4/src/interface/observer_xbraid.h"

namespace ug { namespace XBraidPoroelasticity {

    /**
     * A specialized observer class to write out the result of the BarryMercerProblem as a parallel io gridfunction
     * @tparam TDomain
     * @tparam TAlgebra
     */
    template <typename TDomain, typename TAlgebra>
    class BraidBiotCheck : public XBraidForUG4::IXBraidTimeIntegratorObserver<TDomain, TAlgebra> {
    public:

        //--------------------------------------------------------------------------------------------------------------

        typedef GridFunction<TDomain, TAlgebra> T_GridFunction;
        typedef SmartPtr<T_GridFunction> SP_GridFunction;

        typedef Poroelasticity::BarryMercerProblem<TDomain, TAlgebra> T_Problem;
        typedef SmartPtr<T_Problem> SP_Problem;

        typedef XBraidForUG4::PIOGridFunction<TDomain, TAlgebra> T_PIOGridFunction;
        typedef SmartPtr<T_PIOGridFunction> SP_PIOGridFunction;

        //--------------------------------------------------------------------------------------------------------------

        SP_Problem m_problem;
        int napprox = 16;
        const char * filename = "BarryMercer_2D";

        //--------------------------------------------------------------------------------------------------------------

        BraidBiotCheck() : XBraidForUG4::IXBraidTimeIntegratorObserver<TDomain, TAlgebra>() {};

        ~BraidBiotCheck() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void set_problem(SP_Problem problem) {
            this->m_problem = problem;
        }

        void set_napprox(int napprox) {
            this->napprox = napprox;
        }

        void set_filename(const char * filename) {
            this->filename = filename;
        }

        virtual bool step_process(SP_GridFunction u, int index, double time, double dt) {
            m_problem->m_errData.napprox = this->napprox;

            SP_GridFunction solution = m_problem->compute_solution(u, index, time);

            std::stringstream ss;
            ss << this->filename<<"_"<<index;
            SP_PIOGridFunction pio = SP_PIOGridFunction();
            pio->write(solution,ss.str().c_str());

            return false;
        };

        virtual bool step_process(SP_GridFunction u, int index, double time, double dt, int level, int interation) {
            this->step_process(u,index,time,dt);
            return false;
        };


        bool lua_write(SP_GridFunction u, int index, double time) {
            return step_process(u, index, time, 0);
        }
        //--------------------------------------------------------------------------------------------------------------
    };
}}


#endif //UG_PLUGIN_XBRAIDBIOT_BRAID_BIOT_CONTROL_H
