#ifndef EXTRANODEINFO_H_
#define EXTRANODEINFO_H_

// utility class to allow SearchNodes to store different information
// depending on the heuristic class

namespace daoopt {

class ExtraNodeInfo {
public:
    virtual ~ExtraNodeInfo() = 0;
};

inline ExtraNodeInfo::~ExtraNodeInfo() { }

}  // namespace daoopt

#endif
