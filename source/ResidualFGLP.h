#ifndef _RESIDUALFGLP_H_
#define _RESIDUALFGLP_H_

#include "Problem.h"
#include "Function.h"
#include "mex/indexedHeap.h"
#include <map>

// Implementation of Residual FGLP.

// Class to contain all of the f->i and i->f messages for a variable
// as well as the nullary shift (if applicable)
// also stores the *previous* i->f message
namespace daoopt {

enum Distance { L1, L2, LInf };
class ResidualFGLP;
class FGLPVariableUpdate {
    friend class ResidualFGLP;
    public:
        inline FGLPVariableUpdate() 
            : var_(-1), domain_(-1), 
            prev_var_to_f_(nullptr), var_to_f_(nullptr), 
            prev_nullary_shift_(ELEM_ONE), nullary_shift_(ELEM_ONE),
            owner_(nullptr) {}

        inline FGLPVariableUpdate(int var, int domain) 
            : var_(var), domain_(domain), 
            prev_nullary_shift_(ELEM_ONE), nullary_shift_(ELEM_ONE),
            owner_(nullptr) {
                prev_var_to_f_ = new double[domain];
                var_to_f_ = new double[domain];
                for (int i = 0; i < domain; ++i) {
                    prev_var_to_f_[i] = ELEM_ONE;
                    var_to_f_[i] = ELEM_ONE;
                }
            }

        // Adds or updates the message that the factor sends to the variable
        inline void UpdateFactorMessage(int fid, double *message) {
            if (factor_to_var_.find(fid) != factor_to_var_.end())
                delete factor_to_var_[fid];
            factor_to_var_[fid] = message;
        }


        inline void UpdateVariableMessage(double *message) {
            delete prev_var_to_f_;
            prev_var_to_f_ = var_to_f_;
            var_to_f_ = message;
        }

        inline void ReplaceVariableMessage(double *message) {
            delete var_to_f_;
            var_to_f_ = message;
        }

        inline void UpdateNullaryShift(double shift) {
            prev_nullary_shift_ = nullary_shift_;
            nullary_shift_ = shift;
        }

        inline void ReplaceNullaryShift(double shift) {
            nullary_shift_ = shift;
        }

        // Calculate the residual using the specified distance measure
        // (based on the difference of the variable message and nullary shift)
        inline double UpdateResidual(Distance d=Distance::LInf) {
//            cout << "var: " << var_ << endl;
            residual_ = ELEM_ONE;

            /* DEBUG show previous and current message */
            /*
            for (int i = 0; i < domain_; ++i) {
                cout << prev_var_to_f_[i] << " " << var_to_f_[i] << endl;
            }
            cout << endl;
            */
            switch(d) {
                case Distance::L1:
                    for (int i = 0; i < domain_; ++i) {
                        residual_ += fabs(prev_var_to_f_[i] - var_to_f_[i]); 
                    }
                    residual_ += fabs(prev_nullary_shift_ - nullary_shift_);
                    break;
                case Distance::L2:
                    for (int i = 0; i < domain_; ++i) {
                        residual_ += pow(prev_var_to_f_[i] - var_to_f_[i], 2);
                    }
                    residual_ += pow(prev_nullary_shift_ - nullary_shift_, 2);
                    break;
                case Distance::LInf:
                    residual_ = ELEM_ZERO;
                    for (int i = 0; i < domain_; ++i) {
                        residual_ = max(residual_, fabs(prev_var_to_f_[i] - var_to_f_[i]));
                    }
                    residual_ += fabs(prev_nullary_shift_ - nullary_shift_);

                    /*
                    cout << "r = " << residual_
                        << " (" << mes_r << "+" << shift_r << ")" << endl;
                        */
                    break;
                default:
                    throw std::runtime_error("Invalid distance type");
            }
            return residual_;
        }

        inline double nullary_shift() const { return nullary_shift_; }
        inline double residual() const { return residual_; }

        inline void set_owner(ResidualFGLP *owner) {
            owner_ = owner;
        }

        inline FGLPVariableUpdate *Clone() {
            FGLPVariableUpdate *clone = new FGLPVariableUpdate();
            clone->var_ = var_;
            clone->domain_ = domain_;
            double *table;
            double *old_table;

            // Copy factor to variable messages
            for (const auto &pair : factor_to_var_) {
                table = new double[domain_];
                old_table = pair.second;
                for (int i = 0; i < domain_; ++i) {
                    table[i] = old_table[i];
                }
                clone->factor_to_var_[pair.first] = table;
            }

            // Copy previous variable to factor message
            table = new double[domain_];
            for (int i = 0; i < domain_; ++i) {
                table[i] = prev_var_to_f_[i];
            }
            clone->prev_var_to_f_ = table;

            // Copy current variable to factor message
            table = new double[domain_];
            for (int i = 0; i < domain_; ++i) {
                table[i] = var_to_f_[i];
            }
            clone->var_to_f_ = table;
            
            // Copy nullary shifts
            clone->prev_nullary_shift_ = prev_nullary_shift_;
            clone->nullary_shift_ = nullary_shift_;

            clone->residual_ = residual_;
            return clone;
        }

        inline ~FGLPVariableUpdate() {
            for (auto e : factor_to_var_) delete e.second;
            delete prev_var_to_f_;
            delete var_to_f_;
        }
    private:
        int var_; 
        int domain_;
        std::map<int,double*> factor_to_var_;
        double *prev_var_to_f_;
        double *var_to_f_;
        double prev_nullary_shift_;
        double nullary_shift_;

        double residual_;

        ResidualFGLP *owner_;

};

class ResidualFGLP {
    public:
        static constexpr double DEFAULT_RESIDUAL_TOLERANCE = 1e-8;

        ResidualFGLP(Problem *p, bool useNullaryShift = false);

        ResidualFGLP(ResidualFGLP *parentFGLP, const map<int,val_t> &assignment);
        
        virtual void Run(int maxUpdates, double maxTime, 
                double tolerance=DEFAULT_RESIDUAL_TOLERANCE);

        void GetVarUB(int var, vector<double> &out);

        inline double ub() const { return ub_; }

        inline double global_constant() const { 
            return global_const_factor_->getTable()[0];
        }

        void GetLabelAll(int var, vector<double> &out);
        size_t GetSize() const;

        inline void set_verbose(bool v) { verbose_ = v; }
        inline void set_owns_factors(bool o) { owns_factors_ = o; }

        inline const vector<Function*> &factors() const {
            return factors_;
        }

        inline double runtime() const { return runtime_; }
        inline double runiters() const { return runiters_; }

        inline const set<int> &vars_updated() const {
            return vars_updated_;
        }


        inline virtual ~ResidualFGLP() {
            if (owns_factors_)
                for (Function *f : factors_) delete f;

            for (FGLPVariableUpdate * f_update : message_updates_)
                if (f_update->owner_ == this) delete f_update;
        }

    private:
        // By default, assumes the update was applied and updates the 
        // messages stored to rotate out the previous message
        // In replace mode, this is used for case where the update was not
        // actually applied, but required due to another update that was 
        // applied
        virtual void ComputeVariableMessage(int v, bool replace_mode=false);

        void Condition(const vector<Function*> &fns, 
                const map<int,val_t> &assn);

        double UpdateUB();

        double *MaxMarginal(Function *f, int v);

        void Reparameterize(Function *f, double *maxMarg, double *avgMaxMarg,
                int v);

        Problem *problem_;

        bool owns_factors_;

        vector<Function*> factors_;

        vector<vector<Function*>> factors_by_variable_;

        Function *global_const_factor_;

        set<int> vars_to_update_;

        double ub_;

        bool use_nullary_shift_;

        bool verbose_;

        double runtime_;

        int runiters_;

        set<int> vars_updated_;

        // Priority queue to schedule updates
        vector<FGLPVariableUpdate*> message_updates_;
        mex::indexedHeap message_queue_; 

        vector<int> distance_;

                        
};

}  // namespace daoopt

#endif
