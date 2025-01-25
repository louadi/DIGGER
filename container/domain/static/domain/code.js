var cy = cytoscape({
  container: document.getElementById('cy'),

  boxSelectionEnabled: false,
  autounselectify: true,

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'content': 'data(name)',
        'text-valign': 'center',
        'color': 'white',
        'text-outline-width': 2,
        'text-outline-color': '#888',
        'background-color': '#888'
      })
    .selector(':selected')
      .css({
        'background-color': 'black',
        'line-color': 'black',
        'target-arrow-color': 'black',
        'source-arrow-color': 'black',
        'text-outline-color': 'black'
      }),

      'elements': {'edges': [{'data': {'source': '3670/PF00046',
           'target': '3670/PF00046'}},
         {'data': {'source': '3670/PF00046', 'target': '89884/PF00046'}},
         {'data': {'source': '3670/PF00046', 'target': '64843/PF00046'}},
         {'data': {'source': '64843/PF00046', 'target': '89884/PF00046'}}],
        'nodes': [{'data': {'id': '3670/PF00046',
           'name': '3670/PF00046',
           'value': '3670/PF00046'}},
         {'data': {'id': '64843/PF00046',
           'name': '64843/PF00046',
           'value': '64843/PF00046'}},
         {'data': {'id': '89884/PF00046',
           'name': '89884/PF00046',
           'value': '89884/PF00046'}}]},
      
        layout: {
          name: 'grid',
          padding: 10
        }
      });

cy.on('tap', 'node', function(){
  try { // your browser may block popups
    window.open( this.data('href') );
  } catch(e){ // fall back on url change
    window.location.href = this.data('href');
  }
});