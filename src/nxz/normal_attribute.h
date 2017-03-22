#ifndef NX_NORMAL_ATTRIBUTE_H
#define NX_NORMAL_ATTRIBUTE_H

#include "attribute.h"
#include "point.h"

namespace nx {

class NormalAttr: public GenericAttr<int> {
public:
	NormalAttr(): GenericAttr<int>(2) {}
};

} //namespace
#endif // NX_NORMAL_ATTRIBUTE_H
