#ifndef PLOTTINGXMLDOCSINGLE_HXX_SEEN
#define PLOTTINGXMLDOCSINGLE_HXX_SEEN

#include "TXMLEngine.h"

class EXMLDocAlreadyLoaded : public std::exception {};
class EXMLDocFailedLoad : public std::exception {};
class EXMLDocGetFailed : public std::exception {};
class EXMLDocNullSingleton : public std::exception {};


class XMLDoc {

private:
  XMLDoc(){}
  ~XMLDoc();
  std::map<std::string, TXMLEngine*> docs;

public:
  static void LoadDoc(std::string docID, std::string FileName);
  static TXMLEngine& Get(std::string docID);
  static TXMLEngine& Destroy(std::string docID);

};

class EXMLUtilsFailedToFindElement : public std::exception {};

namespace XMLUtils {
  XMLDocPointer_t GetNamedChildElementOfDocumentRoot(TXMLEngine &xmlengine,
    std::string const & ElementName);

  std::string GetAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName);

  double GetDoubleAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName, bool rethrow=false, double def=0xdeadbeef);

  int GetIntAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName, bool rethrow=false, int def=0);

  bool GetBoolAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName, bool& vaild);


  std::vector<XMLDocPointer_t> GetNamedChildNodes(TXMLEngine &xmlengine,
    XMLDocPointer_t node, std::string const &Name);

  XMLDocPointer_t GetFirstNamedChildNode(TXMLEngine &xmlengine,
    XMLDocPointer_t node, std::string const &Name);

  std::string GetNodeContent(TXMLEngine &xmlengine, XMLDocPointer_t node);

  std::string GetChildNodeContent(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &ChildNodeName);
}

#endif
