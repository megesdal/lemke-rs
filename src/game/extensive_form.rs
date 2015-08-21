
struct Move {
    id: usize,
    label: String,
}

impl Move {
    fn new(id: usize, label: String) -> Move {
        Move { id: id, label: label }
    }
}

struct Node<'a> {
    id: usize,
    sibling: Option<&'a Node<'a>>,
    next_in_iset: Option<&'a Node<'a>>,
}

impl<'a> Node<'a> {
    fn new(id: usize) -> Node<'a> {
        Node {
            id: id,
            sibling: None,
            next_in_iset: None,
        }
    }
}

struct Outcome<'a> {
    id: usize,
    node: &'a Node<'a>,
}

impl<'a> Outcome<'a> {
    fn new(id: usize, node: &'a Node<'a>) -> Outcome<'a> {
        Outcome { id: id, node: node }
    }
}

#[derive(Hash, Eq, PartialEq, Debug)]
struct Player {
    id: usize,
    pub name: String
}

impl Player {

	fn new(id: usize, name: String) -> Player {
        Player { id: id + 1, name: name }
    }

    fn chance() -> Player {
        Player { id: 0, name: "!".to_string() }
    }

    fn is_chance(&self) -> bool {
        return self.id == 0;
    }
}

impl ToString for Player {
	fn to_string(&self) -> String {
		self.name.clone()
	}
}

struct InformationSet<'a> {
    id: usize,
    pub name: String,
    player: &'a Player,
    nodes: Vec<&'a Node<'a>>,
    moves: Vec<&'a Move>,
}

impl<'a> InformationSet<'a> {

    fn new(id: usize, pl: &'a Player) -> InformationSet<'a> {
        InformationSet {
            id: id,
            name: id.to_string(),
            nodes: Vec::new(),
            moves: Vec::new(),
            player: pl
        }
    }
    /*fn name(&self) -> String {
        match self.name {
            Some(x) => x,
            None    => self.id.to_string()
        }
    }*/

	/*fn move_count(&self) -> usize {
    	let mut count = 0;
        let mut child = self.first_node;
        loop {
            if child == None {  // can I do this?  Or do I need a match...
                break;
            }
            count += 1;
            child = child.sibling;
        }
    	count
	}*/

	fn insert_node(&'a mut self, node: &'a Node) {
        // TODO: validate node children equal number of moves?
        // Or lookup the child based on the move?
        self.nodes.push(node);
	}
}

impl<'a> ToString for InformationSet<'a> {
    fn to_string(&self) -> String {
        self.name.clone()
    }
}

struct ExtensiveForm<'a> {

    pub players: Vec<Player>,
    pub chance_player: Player,

    //pub root: &Node,
    pub nodes: Vec<Node<'a>>,

    outcomes: Vec<Outcome<'a>>,

    isets: Vec<InformationSet<'a>>,

    moves: Vec<Move>,
}

impl<'a> ExtensiveForm<'a> {

    fn new() -> ExtensiveForm<'a> {
        ExtensiveForm {
            players: Vec::new(),
            chance_player: Player::chance(),
            nodes: vec![Node::new(0)],
            outcomes: Vec::new(),
            isets: Vec::new(),
            moves: Vec::new(),
        }
    }

    fn root_node(&'a self) -> &'a Node {
        &self.nodes[0]
    }

    fn create_node(&'a mut self) -> &'a Node {
        let id = self.nodes.len();
        self.nodes.push(Node::new(id));
        &self.nodes[id]
    }

    // TODO...
    /*fn first_leaf() -> &Node {
    	root.firstLeaf()
    }*/

    // TODO...
    /* public void autoname()
    {
        for (Player pl = _firstPlayer; pl != null; pl = pl.next)    /* name isets of player pl      */
    	{
        	int idx = pl == Player.CHANCE ? 0 : pl == _firstPlayer ? 1 : 2;
    	    int anbase = an2[idx]-an1[idx]+1;

    	    int digits = 1;
    	    for (int max = anbase, n = nisets(pl); max < n; max *= anbase) {
    	        ++digits;
    	    }

    	    int count = 0;
    	    for (Iset h = _root.iset(); h != null; h = h.next())
    	    {
    	    	if (h.player() == pl) {
	                StringBuilder sb = new StringBuilder();
	        	    for (int j = digits - 1, i = count; j >= 0; --j, i /= anbase)
	        		{
	                    char c = (char)(an1[idx] + (i % anbase));
	            		sb.append(c);
	        		}
	                h.setName(sb.toString());
    	    	}
    	    	++count;
    	    }
    	}
    }*/

    fn num_isets(&self, pl: &Player) -> usize {
        self.isets.iter().filter(|h| h.player == pl).fold(0, |acc, _| acc + 1)
    }

	fn create_outcome(&'a mut self, leaf_node: &'a Node) -> &'a Outcome {
        let outcome_id = self.outcomes.len();
		let new_outcome = Outcome::new(outcome_id, leaf_node);
        self.outcomes.push(new_outcome);
        &self.outcomes[outcome_id]
	}

	fn create_information_set(&'a mut self, name: Option<String>, player: &'a Player) -> &'a InformationSet {

        let iset_id = self.isets.len();
		let h = InformationSet::new(iset_id, player);
		/*if let Some(x) = name {
			h.name = x;
		}*/
		/*if (_secondIset == null) {
			_secondIset = h;
		}
		if let Some(x) = self.last_iset {
			last_iset.next = .setNext(h);
		}
		h.setNext(null);
		_lastIset = h;*/
        self.isets.push(h);
		&self.isets[iset_id]
	}

	fn create_player(&'a mut self, player_name: String) -> &'a Player {
		if player_name == self.chance_player.name {
			&self.chance_player
		} else {
    		let player_id = self.players.len();
            self.players.push(Player::new(player_id, player_name));
    		&self.players[player_id]
        }
	}

    // TODO...
	/*fn add_to_iset(&mut self, node: &mut Node, iset: &mut Iset) {

		node.iset(iset);
		iset.insert_node(node);
		if (node == root) {
			// pull iset out of list & make it the front
			for (Iset h = _secondIset; h != null; h = h.next()) {
				if h.next == iset {
					h.next = iset.next;
				}
			}
			if (iset != second_iset) { //avoid the infinite loop
				iset.next = second_iset;
			}
		}
	}*/

	fn create_move(&'a mut self, move_label: String) -> &'a Move {
		let move_id = self.moves.len();
		self.moves.push(Move::new(move_id, move_label));
		&self.moves[move_id]
	}
}
