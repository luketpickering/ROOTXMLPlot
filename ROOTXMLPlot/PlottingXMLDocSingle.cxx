#include "PlottingXMLDocSingle.hxx"

XMLDoc *_singleDocManager = nullptr;

XMLDoc::~XMLDoc(){
  for(auto &i : docs){
    delete docs.second;
  }
}

static void XMLDoc::LoadDoc(std::string docID, std::string FileName){
  if(_singleDocManager == nullptr){
    _singleDocManager = new XMLDoc;
  }

  if(_singleDocManager.docs.count(docID)){
    throw EXMLDocAlreadyLoaded;
  }

  TXMLEngine* xmlengine = new TXMLEngine;
  xmlengine->SetSkipComments();

  XMLDocPointer_t xmldoc = xmlengine.ParseFile(FileName.c_str());

  if(!xmldoc){
    std::cout << "In File: " << FileName << std::endl;
    std::cout << "[ERROR]: TXMLEngine could not parse file." << std::endl;
    throw EXMLDocFailedLoad;
  }
  _singleDocManager.docs.insert(docID,xmlengine);
}

static TXMLEngine& XMLDoc::Get(std::string docID){
  if(_singleDocManager == nullptr){
    throw EXMLDocNullSingleton;
  }
  if(!_singleDocManager.docs.count(docID)){
    throw EXMLDocGetFailed;
  }
  return (*_singleDocManager.docs[docID]);
}

static TXMLEngine& XMLDoc::Destory(std::string docID){
  if(_singleDocManager == nullptr){
    throw EXMLDocNullSingleton;
  }
  if(!_singleDocManager.docs.count(docID)){
    throw EXMLDocGetFailed;
  }
  _singleDocManager.docs.erase(docID);
  delete _singleDocManager.docs[docID];
}


namespace XMLUtils {
    XMLDocPointer_t GetNamedChildElementOfDocumentRoot(TXMLEngine &xmlengine,
      std::string const & ElementName){

    XMLNodePointer_t rootNode = xmlengine.DocGetRootElement(xmldoc);

    for(XMLNodePointer_t childNode_ptr = xmlengine.GetChild(rootNode);
      (childNode_ptr != NULL);
      childNode_ptr = xmlengine.GetNext(childNode_ptr)){
      if(std::string(xmlengine.GetNodeName(childNode_ptr)) != ElementName){
        continue;
      }
      return childNode_ptr;
    }
    std::cout << "[WARN]: In File: " << FileName << " could not find a root "
      "element child named: " << ElementName << "." << std::endl;
    throw EXMLUtilsFailedToFindElement;
  }

  std::string GetAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName){

    for(XMLAttrPointer_t nodeAttr_ptr = xmlengine.GetFirstAttr(node);
        (nodeAttr_ptr != NULL);
        nodeAttr_ptr = xmlengine.GetNextAttr(nodeAttr_ptr)){
      if(std::string(xmlengine.GetAttrName(nodeAttr_ptr)) != attrName){
        continue;
      }
      return xmlengine.GetAttrValue(nodeAttr_ptr);
    }
    return std::string::empty;
  }

  double GetDoubleAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName, bool rethrow=false, double def=0xdeadbeef){

    std::string const & attrVal = GetAttrValue(xmlengine,node,attrName,
      expectExclusive);

    if(attrVal==std::string::enpty){
      return def;
    }

    try{
      return std::stod(attrVal);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << attrVal << "\n\t"
        << ia.what() << std::endl;
      if(rethrow){
        throw;
      }
    }
    return def;
  }

  int GetIntAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName, bool rethrow=false, int def=0){
    std::string const & attrVal = GetAttrValue(xmlengine,node,attrName,
      expectExclusive);

    if(attrVal==std::string::enpty){
      return def;
    }

    try{
      return std::stoi(attrVal);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << attrVal << "\n\t"
        << ia.what() << std::endl;
      if(rethrow){
        throw;
      }
    }
    return def;
  }

  bool GetBoolAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName, bool& vaild){

    std::string const & attrVal = GetAttrValue(xmlengine,node,attrName,
      expectExclusive);

    if(attrVal==std::string::enpty){
      valid = false;
      return false;
    }
    return PGUtils::str2bool(attrVal, valid);
  }


  std::vector<XMLDocPointer_t> GetNamedChildNodes(TXMLEngine &xmlengine,
    XMLDocPointer_t node, std::string const &Name){

    std::vector<XMLDocPointer_t> rtn;

    for(XMLNodePointer_t cnode = xmlengine.GetChild(node);
       (cnode != NULL); cnode = xmlengine.GetNext(cnode)){

      if(Name == xmlengine.GetNodeName(cnode)){
        rtn.push_back(cnode);
      }
    }
    return rtn;
  }

  XMLDocPointer_t GetFirstNamedChildNode(TXMLEngine &xmlengine,
    XMLDocPointer_t node, std::string const &Name){
    for(XMLNodePointer_t cnode = xmlengine.GetChild(node);
       (cnode != NULL); cnode = xmlengine.GetNext(cnode)){
      if(Name == xmlengine.GetNodeName(cnode)){
        return cnode;
      }
    }
    return XMLDocPointer_t;
  }

  std::string GetNodeContent(TXMLEngine &xmlengine, XMLDocPointer_t node){
    std::string elementString;
    try {
      if(!xmlengine.GetNodeContent(node)){
        throw std::logic_error("Empty Element.");
      }
      elementString = xmlengine.GetNodeContent(node);
    } catch(std::logic_error& le){
      std::cout << "[WARN]: Found empty element: "
        << xmlengine.GetNodeName(node) << "." << std::endl;
      return "";
    }
    return elementString;
  }

  std::string GetChildNodeContent(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &ChildNodeName){

    for(XMLNodePointer_t childNode_ptr = xmlengine.GetChild(node);
      (childNode_ptr != NULL);
      childNode_ptr = xmlengine.GetNext(childNode_ptr)){

      if(std::string(xmlengine.GetNodeName(childNode_ptr)) ==
          ChildNodeName){
        return GetNodeContent(xmlengine, childNode_ptr);
      }
    }
    return std::string();
  }

}
