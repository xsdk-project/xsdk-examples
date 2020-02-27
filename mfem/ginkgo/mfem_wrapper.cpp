#include "mfem_wrapper.hpp"

#include <ginkgo/ginkgo.hpp>

void MFEMOperatorWrapper::apply_impl(const gko::LinOp *b, gko::LinOp *x) const
{

    // Cast to MFEMVectorWrapper; only accept this type for this impl
    const MFEMVectorWrapper *mfem_b = gko::as<const MFEMVectorWrapper>(b);
    MFEMVectorWrapper *mfem_x = gko::as<MFEMVectorWrapper>(x);

    this->mfem_oper_->Mult(mfem_b->get_mfem_vec_const_ref(),
                           mfem_x->get_mfem_vec_ref());
}
void MFEMOperatorWrapper::apply_impl(const gko::LinOp *alpha,
                                     const gko::LinOp *b,
                                     const gko::LinOp *beta,
                                     gko::LinOp *x) const
{
    // x = alpha * op (b) + beta * x

    // Cast to MFEMVectorWrapper; only accept this type for this impl
    const MFEMVectorWrapper *mfem_b = gko::as<const MFEMVectorWrapper>(b);
    MFEMVectorWrapper *mfem_x = gko::as<MFEMVectorWrapper>(x);

    // Check that alpha and beta are Dense<double> of size (1,1):
    if (alpha->get_size()[0] > 1 || alpha->get_size()[1] > 1) {
        throw gko::BadDimension(
            __FILE__, __LINE__, __func__, "alpha", alpha->get_size()[0],
            alpha->get_size()[1],
            "Expected an object of size [1 x 1] for scaling "
            " in this operator's apply_impl");
    }
    if (beta->get_size()[0] > 1 || beta->get_size()[1] > 1) {
        throw gko::BadDimension(
            __FILE__, __LINE__, __func__, "beta", beta->get_size()[0],
            beta->get_size()[1],
            "Expected an object of size [1 x 1] for scaling "
            " in this operator's apply_impl");
    }

    double alpha_f;
    double beta_f;

    if (alpha->get_executor() == alpha->get_executor()->get_master()) {
          // Access value directly
          alpha_f = gko::as<gko::matrix::Dense<double>>(alpha)->at(0, 0);
    } else {
        // Copy from device to host
        this->get_executor()->get_master().get()->copy_from(this->get_executor().get(),
              1, gko::as<gko::matrix::Dense<double>>(alpha)->get_const_values(),
              &alpha_f); 
    }

    if (beta->get_executor() == beta->get_executor()->get_master()) {
          // Access value directly
          beta_f = gko::as<gko::matrix::Dense<double>>(beta)->at(0, 0);
    } else {
        // Copy from device to host
        this->get_executor()->get_master().get()->copy_from(this->get_executor().get(),
              1, gko::as<gko::matrix::Dense<double>>(beta)->get_const_values(),
              &beta_f); 
    }

    // Scale x by beta
    mfem_x->get_mfem_vec_ref() *= beta_f;

    // Multiply operator with b and store in tmp
    mfem::Vector mfem_tmp =
        mfem::Vector(mfem_x->get_size()[0],
                     mfem_x->get_mfem_vec_ref().GetMemory().GetMemoryType());

    // Set UseDevice flag to match mfem_x (not automatically done through
    // MemoryType)
    mfem_tmp.UseDevice(mfem_x->get_mfem_vec_ref().UseDevice());

    this->mfem_oper_->Mult(mfem_b->get_mfem_vec_const_ref(), mfem_tmp);

    // Scale tmp by alpha and add
    mfem_x->get_mfem_vec_ref().Add(alpha_f, mfem_tmp);

    mfem_tmp.Destroy();
}
