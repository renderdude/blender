from bs4 import BeautifulSoup
import os, sys

def generate_mtlx(data):
    bs_data = BeautifulSoup(data, "xml")

    node_defs = bs_data.find_all('nodedef')
    for nd in node_defs:
        nd_attrs = dict(nd.attrs)
        out_type = nd.find('output')
        out_attrs = dict(out_type.attrs)
        node_tag = nd_attrs['node']
        if len(node_defs) > 1:
            node_name = node_tag + out_attrs['type']
        else:
            node_name = node_tag
            
        str = '  <!-- ' + node_name + ' -->'
        print(str)
        str = '  <' + node_tag + ' name="' + node_name + '" type="' + out_attrs['type'] + '">'
        print(str)

        params = nd.find_all('input')
        for p in params:
            param_dict = dict(p.attrs)
            str = '    <input name="' + param_dict['name'] + '"'
            str += ' type="' + param_dict['type'] + '"'
            if 'value' in param_dict:
                str += ' value="' + param_dict['value'] + '"'
            str += '/>'
            print(str)

        str = '  </' + node_tag + '>'
        print(str)
        str = '  <surface name="' + node_name + 'Surface" type="surfaceshader">'
        print(str)
        str = '    <input name="' + out_attrs['type'].lower() + '" type="' + out_attrs['type'] + '" nodename="' + node_name + '" />'
        print(str)
        str = '  </surface>'
        print(str)
        str = '  <surfacematerial name="' + node_name + 'Shader" type="material">'
        print(str)
        str = '    <input name="surfaceshader" type="surfaceshader" nodename="' + node_name +'Surface" />'
        print(str)
        str = '  </surfacematerial>\n'
        print(str)

def main():
    if len(sys.argv) != 2:
        print('Need path to MaterialX lama node definitions')
        sys.exit(-1)
    
    print('<?xml version="1.0"?>')
    print('<materialx version="1.38" colorspace="acescg">\n')

    for root, dirs, files in os.walk(sys.argv[1]):
        for name in files:
            if name.endswith('mtlx'):
                with open(os.path.join(root, name), 'r') as f:
                    data = f.read()

                generate_mtlx(data)
    print('</materialx>')

if __name__ == "__main__":
    main()
