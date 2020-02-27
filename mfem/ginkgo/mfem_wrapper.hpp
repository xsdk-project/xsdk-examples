#include "mfem.hpp"

#include <ginkgo/ginkgo.hpp>

template <typename T>
class mfem_destroy {
public:
    using pointer = T *;

    // Destroys an MFEM object.  Requires object to have a Destroy() method.
    void operator()(pointer ptr) const noexcept { ptr->Destroy(); }
};

// This class wraps an MFEM vector object for Ginkgo's use.  It 
// is derived from Ginkgo's Dense class, meaning we can perform
// standard Ginkgo operations for Dense vectors, while also passing
// through the MFEM Vector object to an MFEM Operator.  The Ginkgo
// data is set to non-owning views of the MFEM Vector's data, so 
// no copying is involved and changes made by one library will be 
// seen by the other.
class MFEMVectorWrapper : public gko::matrix::Dense<double> {
public:
    MFEMVectorWrapper(std::shared_ptr<const gko::Executor> exec,
                      gko::size_type size, mfem::Vector *mfem_vec,
                      bool on_device = true, bool ownership = false)
        : gko::matrix::Dense<double>(
              exec, gko::dim<2>{size, 1},
              gko::Array<double>::view(exec, size,
                                       mfem_vec->ReadWrite(on_device)),
              1)
    {
        // This controls whether or not we want Ginkgo to own its MFEM Vector.
        // Normally, when we are wrapping an MFEM Vector created outside 
        // Ginkgo, we do not want ownership to be true. However, Ginkgo
        // creates its own temporary vectors as part of its solvers, and 
        // these will be owned (and deleted) by Ginkgo. 
        if (ownership) {
            using deleter = mfem_destroy<mfem::Vector>;
            mfem_vec_ = std::unique_ptr<mfem::Vector,
                                        std::function<void(mfem::Vector *)>>(
                mfem_vec, deleter{});
        } else {
            using deleter = gko::null_deleter<mfem::Vector>;
            mfem_vec_ = std::unique_ptr<mfem::Vector,
                                        std::function<void(mfem::Vector *)>>(
                mfem_vec, deleter{});
        }
    }

    // The create function returns a unique_ptr to an MFEMVectorWrapper object.
    // This is the function that should be called from the program code.
    static std::unique_ptr<MFEMVectorWrapper> create(
        std::shared_ptr<const gko::Executor> exec, gko::size_type size,
        mfem::Vector *mfem_vec, bool on_device = true, bool ownership = false)
    {
        return std::unique_ptr<MFEMVectorWrapper>(
            new MFEMVectorWrapper(exec, size, mfem_vec, on_device, ownership));
    }

    // Return reference to MFEM Vector object
    mfem::Vector &get_mfem_vec_ref() const { return *(this->mfem_vec_.get()); }
    // Return const reference to MFEM Vector object
    const mfem::Vector &get_mfem_vec_const_ref() const
    {
        return const_cast<const mfem::Vector &>(*(this->mfem_vec_.get()));
    }

    // Override base Dense class implementation for creating new vectors
    // with same executor and size as self
    virtual std::unique_ptr<gko::matrix::Dense<double>>
    create_with_same_config() const override
    {
        mfem::Vector *mfem_vec = new mfem::Vector(
            this->get_size()[0],
            this->mfem_vec_.get()->GetMemory().GetMemoryType());

        mfem_vec->UseDevice(this->mfem_vec_.get()->UseDevice());

        // If this function is called, Ginkgo is creating this 
        // object and should control the memory, so ownership is 
        // set to true
        return MFEMVectorWrapper::create(
            this->get_executor(), this->get_size()[0], mfem_vec,
            this->mfem_vec_.get()->UseDevice(), true);
    }

private:
    std::unique_ptr<mfem::Vector, std::function<void(mfem::Vector *)>>
        mfem_vec_;
};

// This class wraps an MFEM Operator for Ginkgo, to make its Mult()
// function available to Ginkgo, provided the input and output vectors
// are of MFEMVectorWrapper type.
//
// It is also an example of using `EnableLinOp` to create custom
// linear operators in Ginkgo. 
class MFEMOperatorWrapper
    : public gko::EnableLinOp<MFEMOperatorWrapper>,
      public gko::EnableCreateMethod<MFEMOperatorWrapper> {
public:
    MFEMOperatorWrapper(std::shared_ptr<const gko::Executor> exec,
                        gko::size_type size = 0,
                        mfem::OperatorHandle oper = mfem::OperatorHandle())
        : gko::EnableLinOp<MFEMOperatorWrapper>(exec, gko::dim<2>{size}),
          gko::EnableCreateMethod<MFEMOperatorWrapper>()
    {
        this->mfem_oper_ = oper;
    }
 
protected:
    // These are the two functions that custom Ginkgo LinOps must implement
    void apply_impl(const gko::LinOp *b, gko::LinOp *x) const override;
    void apply_impl(const gko::LinOp *alpha, const gko::LinOp *b,
                    const gko::LinOp *beta, gko::LinOp *x) const override;

private:
    mfem::OperatorHandle mfem_oper_;
};
