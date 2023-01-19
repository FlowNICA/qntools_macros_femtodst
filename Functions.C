// Function to read either Qn::DataContainerStatCollect or Qn::DataContainerStatCalculate from the file into dc_calc
Bool_t GetDCStatCalculate(TFile *const& file, std::string name, Qn::DataContainerStatCalculate& dc_calc)
{
    Qn::DataContainerStatCollect *tmp_dc_collect{nullptr};
    Qn::DataContainerStatCalculate *tmp_dc_calculate{nullptr};

    file->GetObject(name.c_str(), tmp_dc_calculate);
    if (!tmp_dc_calculate)
    {
        file->GetObject(name.c_str(), tmp_dc_collect);
        if (!tmp_dc_collect)
        {
            std::cerr << "Cannot get the object " << name << "!" << std::endl;
            return false;
        }
        dc_calc = (Qn::DataContainerStatCalculate) *tmp_dc_collect;
    }
    else
    {
        dc_calc = (Qn::DataContainerStatCalculate) *tmp_dc_calculate;
    }

    return true;
}

// Returns mean value from dc1 and dc2
Qn::DataContainerStatCalculate Merge(Qn::DataContainerStatCalculate dc1, Qn::DataContainerStatCalculate dc2)
{
    Qn::DataContainerStatCalculate result = dc1;
    if (dc1.size() != dc2.size())
    {
        std::cerr << "Merge: dc1 and dc2 has different sizes!" << std::endl;
        return result;
    }

    for (int i=0; i<result.size(); i++)
    {
        result.At(i) = Qn::Merge(dc1.At(i), dc2.At(i));
    }

    return result;
}